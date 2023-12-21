#include <iostream>
#include <limits>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

constexpr int numnodes    = 13900;
constexpr int numelements = 27307;

//constexpr int numnodes    = 2366;
//constexpr int numelements = 4185;

constexpr double p  = 10;
constexpr double E  = 100;
constexpr double nu = 0.25;

int main(){
    std::ifstream in;
    std::vector<Eigen::Vector2d>          nodes(numnodes);
    std::vector<Eigen::Vector2d>          pOnNodes(numnodes, {0,0});
    std::vector<std::size_t>              fixNodesX;
    std::vector<std::size_t>              fixNodesY;
    std::vector<std::vector<std::size_t>> elements(numelements);

    //in.open("../Data_0.5.k");
    in.open("../Data.k");
    std::cout << "Opened file" << std::endl;
    for(int i = 0; i < 2; i++){
        in.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }

    double ymax = 0;

    for(int i = 0; i < numnodes; i++){
        unsigned int nid;
        unsigned int tmpskip;
        double x,y,z;
        in >> nid >> x >> y >> z >> tmpskip >> tmpskip;
        Eigen::Vector2d v(x,y);

        ymax = std::max(ymax,y);

        nodes[nid-1] = v;
    }
    std::cout << "ymax: " << ymax << std::endl;
    for(int i = 0; i < numnodes; i++){
        if     (std::abs(nodes[i](0)) < std::numeric_limits<double>::min()) fixNodesX.push_back(i);
        else if(std::abs(nodes[i](1)) < std::numeric_limits<double>::min()) fixNodesY.push_back(i);

        if(std::abs(nodes[i](1) - ymax) < std::numeric_limits<double>::min()) {
            pOnNodes[i](1) = p;
        }
    }

    std::cout << fixNodesX.size() + fixNodesY.size() << std::endl;

    for(int i = 0; i < 5; i++){
        in.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }

    for(int i = 0; i < numelements; i++){
        unsigned int eid;
        unsigned int tmpskip;
        unsigned int n1,n2,n3;
        in >> eid >> tmpskip >> n1 >> n2 >> n3 >> tmpskip;
        elements[eid-1].push_back(n1-1);
        elements[eid-1].push_back(n2-1);
        elements[eid-1].push_back(n3-1);
    }

    in.close();

    std::cout << "Loaded mesh and boundary conditions" << std::endl;

    std::vector<Eigen::Vector3d> N(3 * numelements);

    for(int i = 0; i < numelements; i++){
        Eigen::Matrix3d A;
        Eigen::Vector3d b;
        A << 1, nodes[elements[i][0]](0), nodes[elements[i][0]](1), 
             1, nodes[elements[i][1]](0), nodes[elements[i][1]](1), 
             1, nodes[elements[i][2]](0), nodes[elements[i][2]](1);
        b << 1,0,0;
        N[3 * i] = A.colPivHouseholderQr().solve(b);
        b << 0,1,0;
        N[3 * i + 1] = A.colPivHouseholderQr().solve(b);
        b << 0,0,1;
        N[3 * i + 2] = A.colPivHouseholderQr().solve(b);
    }

    Eigen::SparseMatrix<double> Kglob(2*numnodes,2*numnodes);

    for(int i = 0; i < numelements; i++){
        Eigen::Vector3d v1(nodes[elements[i][0]][0],nodes[elements[i][0]][1],0);
        Eigen::Vector3d v2(nodes[elements[i][1]][0],nodes[elements[i][1]][1],0);
        Eigen::Vector3d v3(nodes[elements[i][2]][0],nodes[elements[i][2]][1],0);

        Eigen::Vector3d vv1 = v3 - v1;
        Eigen::Vector3d vv2 = v2 - v1;

        double S = 0.5 * vv1.cross(vv2).norm();

        Eigen::Matrix<double,3,6> B;
        for(int j = 0; j < 3; j++){
            B(0,2*j  ) = N[3 * i + j](1);
            B(0,2*j+1) = 0;
            B(1,2*j  ) = 0;
            B(1,2*j+1) = N[3 * i + j](2);
            B(2,2*j  ) = N[3 * i + j](2);
            B(2,2*j+1) = N[3 * i + j](1);
        }
        /*
        B(0,2) = v3[1] - v1[1];
        B(0,3) = 0;
        B(1,2) = 0;
        B(1,3) = v1[0] - v3[0];
        B(2,2) = v1[0] - v3[0];
        B(2,3) = v3[1] - v1[1];

        B(0,4) = v1[1] - v2[1];
        B(0,5) = 0;
        B(1,4) = 0;
        B(1,5) = v2[0] - v1[0];
        B(2,4) = v2[0] - v1[0];
        B(2,5) = v1[1] - v2[1];

        B = B / (2 * S);
        */
        Eigen::Matrix3d D;
        D << 1, nu / (1 - nu), 0,
        nu / (1 - nu), 1, 0,
        0, 0, (1 - 2 * nu) / (2 * (1 - nu));

        D = D * (E * (1 - nu) / ((1 + nu)*(1 - 2 * nu)));

        Eigen::Matrix<double, 6, 6> Kloc;
        Kloc = B.transpose() * D * B * S;

        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                Kglob.coeffRef(2 * elements[i][j]    , 2 * elements[i][k])     += Kloc(2*j    , 2*k    );
                Kglob.coeffRef(2 * elements[i][j] + 1, 2 * elements[i][k])     += Kloc(2*j + 1, 2*k    );
                Kglob.coeffRef(2 * elements[i][j]    , 2 * elements[i][k] + 1) += Kloc(2*j    , 2*k + 1);
                Kglob.coeffRef(2 * elements[i][j] + 1, 2 * elements[i][k] + 1) += Kloc(2*j + 1, 2*k + 1);
            }
        }
    }
#ifdef CHECK_K
    std::cout << "Checkin Kglob" << std::endl;

    for(int i = 0; i < 2 * numnodes; i++){
        if(abs(Kglob.row(i).sum()) >= 1e-10){
            std::cout << "ERROR ROW NOT NIL!: " << i << " " << Kglob.row(i).sum() << std::endl;
        }
        if(abs(Kglob.col(i).sum()) >= 1e-10){
            std::cout << "ERROR COL NOT NIL!: " << i << " " << Kglob.col(i).sum() << std::endl;
        }
        if(Kglob.coeff(i,i) < 0){
            std::cout << "ERROR DIAG NOT POSITIVE!: " << i << " " << Kglob.coeff(i,i) << std::endl;
        }
    }
    for(int i = 0; i < 2 * numnodes; i++){
        for(int j = 0; j < i; j++){
            if(abs(Kglob.coeff(i, j) - Kglob.coeff(j, i)) >= 1e-10){
                std::cout << "ERROR NON SYMMETRIC!: " << i << " " << j << " " << Kglob.coeff(i,j) << " " << Kglob.coeff(j,i) << std::endl;
            }
        }
    }
#endif

    Eigen::SparseVector<double> Fglob(2*numnodes);
    
    int counter = 0;

    for(int i = 0; i < numelements; i++){
        bool n1 = (pOnNodes[elements[i][0]] != Eigen::Vector2d(0,0));
        bool n2 = (pOnNodes[elements[i][1]] != Eigen::Vector2d(0,0));
        bool n3 = (pOnNodes[elements[i][2]] != Eigen::Vector2d(0,0));
        if(n1 + n2 + n3 < 2) continue;
        if(n1 + n2 + n3 == 3){
            std::cout << "Something weird is going on with P" << std::endl;
        }
        Eigen::Vector<double, 6> Floc;
        Floc(0) = pOnNodes[elements[i][0]][0];
        Floc(1) = pOnNodes[elements[i][0]][1];
        Floc(2) = pOnNodes[elements[i][1]][0]; 
        Floc(3) = pOnNodes[elements[i][1]][1]; 
        Floc(4) = pOnNodes[elements[i][2]][0]; 
        Floc(5) = pOnNodes[elements[i][2]][1];
        double xi1;
        double xi;
        if(n1 && n2){
            xi = nodes[elements[i][0]](0);
            xi1 = nodes[elements[i][1]](0);
        }else if(n1 && n3){
            xi = nodes[elements[i][0]](0);
            xi1 = nodes[elements[i][2]](0);
        }else if(n2 && n3){
            xi = nodes[elements[i][1]](0);
            xi1 = nodes[elements[i][2]](0);
        }

        for(int j = 0; j < 3; j++){
            Floc(2*j+1) *= (ymax * N[3 * i + j](2) + N[3 * i + j](0)) * std::abs(xi1 - xi) +
            N[3 * i + j](1) * 0.5 * std::abs(xi1 - xi) * (xi1 + xi);
        }
        
        Fglob.coeffRef(2*elements[i][0]    ) += Floc(0);
        Fglob.coeffRef(2*elements[i][0] + 1) += Floc(1);
        Fglob.coeffRef(2*elements[i][1]    ) += Floc(2);
        Fglob.coeffRef(2*elements[i][1] + 1) += Floc(3);
        Fglob.coeffRef(2*elements[i][2]    ) += Floc(4);
        Fglob.coeffRef(2*elements[i][2] + 1) += Floc(5);
    }

#ifdef OUT_F
    for(int i = 0; i < numnodes; i++){
        std::cout << Fglob(2*i) << " " << Fglob(2*i + 1) << std::endl;
    }
#endif

    std::cout << "Adding boundary conditions: " << std::endl; 

    for(int i = 0; i < fixNodesX.size(); i++){
        Fglob.coeffRef(2 * fixNodesX[i]) = 0;
        for(int j = 0; j < 2 * numnodes; j++){
            if(Kglob.coeff(2*fixNodesX[i], j) != 0){
                Kglob.coeffRef(2*fixNodesX[i],j             ) = 0;
                Kglob.coeffRef(j,             2*fixNodesX[i]) = 0;
            }
        }
        Kglob.coeffRef(2*fixNodesX[i],2*fixNodesX[i]) = 1;
    }

    for(int i = 0; i < fixNodesY.size(); i++){
        Fglob.coeffRef(2 * fixNodesY[i] + 1) = 0;
        for(int j = 0; j < 2 * numnodes; j++){
            if(Kglob.coeff(2*fixNodesY[i] + 1, j) != 0){
                Kglob.coeffRef(2*fixNodesY[i] + 1,j                 ) = 0;
                Kglob.coeffRef(j,                 2*fixNodesY[i] + 1) = 0;
            }
        }
        Kglob.coeffRef(2*fixNodesY[i] + 1,2*fixNodesY[i] + 1) = 1;
    }
    
    std::cout << "Solver started" << std::endl;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(Kglob);
    Eigen::VectorXd U = ldlt.solve(Fglob);

    std::cout << "Found U." << std::endl;

#ifdef OUT_SIGMA
    std::vector<Eigen::Vector3d> sigma(numelements);

    for(int i = 0; i < numelements; i++){
        Eigen::Matrix<double,3,6> B;
        for(int j = 0; j < 3; j++){
            B(0,2*j  ) = N[3 * i + j](1);
            B(0,2*j+1) = 0;
            B(1,2*j  ) = 0;
            B(1,2*j+1) = N[3 * i + j](2);
            B(2,2*j  ) = N[3 * i + j](2);
            B(2,2*j+1) = N[3 * i + j](1);
        }

        Eigen::Matrix3d D;
        D << 1,             nu / (1 - nu), 0,
             nu / (1 - nu), 1,             0,
             0,             0,             (1 - 2 * nu) / (2 * (1 - nu));

        D = D * (E * (1 - nu) / ((1 + nu)*(1 - 2 * nu)));

        Eigen::Vector<double,6> Uloc;
        Uloc(0) = U(2 * elements[i][0]);
        Uloc(1) = U(2 * elements[i][0]+1);
        Uloc(2) = U(2 * elements[i][1]);
        Uloc(3) = U(2 * elements[i][1]+1);
        Uloc(4) = U(2 * elements[i][2]);
        Uloc(5) = U(2 * elements[i][2]+1);

        Eigen::Vector3d sigmaloc;

        sigmaloc = D * B * Uloc;

        sigma[i](0) = sigmaloc(0);
        sigma[i](1) = sigmaloc(1);
        sigma[i](2) = sigmaloc(2);
    }

    //for(int i = 0; i < numelements; i++){
    //    if(sigma(3 * i + 2) >= 1e-3){
    //        std::cout << sigma(3 * i + 2) << std::endl;
    //    }
    //}

    //for(int i = 0; i < numelements; i++){
    //    std::cout << sigma(3 * i) << " " << sigma(3 * i + 2) << std::endl;
    //}

    Eigen::SparseMatrix<double> Cglob(numnodes, numnodes);
    Eigen::VectorXd Rglobxx(numnodes);
    Eigen::VectorXd Rglobyy(numnodes);
    Eigen::VectorXd Rglobxy(numnodes);
    Rglobxx.setZero();
    Rglobyy.setZero();
    Rglobxy.setZero();

    for(int i = 0; i < numelements; i++){
        Eigen::Matrix3d Cloc;
        Eigen::Vector3d Rloc;
        Eigen::Vector2d a = nodes[elements[i][1]] - nodes[elements[i][0]];
        Eigen::Vector2d b = nodes[elements[i][2]] - nodes[elements[i][0]];

        //Eigen::Matrix4d A;
        //A << a(0), a(1), 0,    0,
        //     0,    0,    a(0), a(1),
        //     b(0), b(1), 0,    0,
        //     0,     0,    b(0), b(1);
        //Eigen::Vector4d ans(1,0,0,1);
        //Eigen::Vector4d perexMat = A.colPivHouseholderQr().solve(ans);
        //Eigen::Matrix2d tmp;
        //tmp << perexMat(0), perexMat(1), perexMat(2), perexMat(3);
        //Eigen::Matrix2d jacob = tmp.inverse();
        Eigen::Matrix2d jacob; 
        jacob << a(0), b(0), a(1), b(1);

        Eigen::Vector3d v1(a(0),a(1),0);
        Eigen::Vector3d v2(b(0),b(1),0);

        double s = 0.5 * v1.cross(v2).norm();

        double s1 = jacob.determinant() * 0.5;

        if(std::abs(s - s1) >= 1e-14){
            std::cout << "AREA CHECK ERROR: " << s << " " << s1 << " " << s - s1 << std::endl;
        }

        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                double a1 = N[3 * i + j](0);
                double b1 = N[3 * i + j](1);
                double c1 = N[3 * i + j](2);

                a1 = a1 + b1 * nodes[elements[i][0]](0) + c1 * nodes[elements[i][0]](1);
                b1 = b1 * jacob(0,0) + c1 * jacob(1,0);
                c1 = b1 * jacob(0,1) + c1 * jacob(1,1);

                double a2 = N[3 * i + k](0);
                double b2 = N[3 * i + k](1);
                double c2 = N[3 * i + k](2);

                a2 = a2 + b2 * nodes[elements[i][0]](0) + c2 * nodes[elements[i][0]](1);
                b2 = b2 * jacob(0,0) + c2 * jacob(1,0);
                c2 = b2 * jacob(0,1) + c2 * jacob(1,1);

                Cloc(j,k) = jacob.determinant() / 24.0 * (4 * a1 * (3 * a2 + b2 + c2) + 4 * a2 * (b1 + c1) + 2 * b1 * b2 + b1 * c2 + b2 * c1 + 2 * c1 * c2);
            }
        }
        for(int j = 0; j < 3; j++){
            double a1 = N[3 * i + j](0);
            double b1 = N[3 * i + j](1);
            double c1 = N[3 * i + j](2);

            a1 = a1 + b1 * nodes[elements[i][0]](0) + c1 * nodes[elements[i][0]](1);
            b1 = b1 * jacob(0,0) + c1 * jacob(1,0);
            c1 = b1 * jacob(0,1) + c1 * jacob(1,1);

            Rloc(j) = jacob.determinant() * 1.0 / 6.0 * (3 * a1 + b1 + c1);
        }
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                Cglob.coeffRef(elements[i][j],elements[i][k]) += Cloc(j,k);
            }
        }
        
        for(int j = 0; j < 3; j++){
            Rglobxx(elements[i][j]) += sigma[i](0) * Rloc(j);
            Rglobyy(elements[i][j]) += sigma[i](1) * Rloc(j);
            Rglobxy(elements[i][j]) += sigma[i](2) * Rloc(j);
        }
    }

    ldlt.compute(Cglob);

    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt2(Cglob);
    Eigen::VectorXd sigmanewxx = ldlt.solve(Rglobxx);
    Eigen::VectorXd sigmanewyy = ldlt.solve(Rglobyy);
    Eigen::VectorXd sigmanewxy = ldlt.solve(Rglobxy);

    double maxeps = 0;

    for(int i = 0; i < numnodes; i++){
        if(nodes[i](1) == 0){
            double analyt = 0.5 * p * (2 + std::pow(1.0 / (nodes[i](0)),2) + 3 * std::pow(1.0 / (nodes[i](0)),4));
            maxeps = std::max(maxeps, std::abs((sigmanewyy(i) - analyt) / analyt));
            /* if(std::abs((sigmanewyy(i) - analyt) / analyt) >= 0.1){
                std::cout << sigmanewyy(i) << " " << analyt << std::endl;
            } */
        }
    }
    
    std::cout << "max error on yy is " << maxeps*100.0 << "%" << std::endl;

    maxeps = 0;

    for(int i = 0; i < numnodes; i++){
        if(nodes[i](1) == 0){
            double analyt = 1.5 * p * (std::pow(1.0 / (nodes[i](0)),2) - std::pow(1.0 / (nodes[i](0)),4));
            if(analyt == 0){
                maxeps = std::max(maxeps, std::abs(sigmanewxx(i)));
                /* if(std::abs(sigmanewxx(i)) >= 0.1){
                    std::cout << sigmanewxx(i) << " " << 0 << std::endl;
                } */
            }else{
                maxeps = std::max(maxeps, std::abs((sigmanewxx(i) - analyt) / analyt));
                /* if(std::abs((sigmanewxx(i) - analyt) / analyt) >= 0.1){
                    std::cout << nodes[i](0) << " " << sigmanewxx(i) << " " << analyt << std::endl;
                } */
            }
        }
    }
    
    std::cout << "max error on xx is " << maxeps*100.0 << "%" << std::endl;

    maxeps = 0;

    for(int i = 0; i < numnodes; i++){
        if(nodes[i](1) == 0){
            maxeps = std::max(maxeps, std::abs(sigmanewxy(i)));
            /* if(std::abs(sigmanewxy(i)) >= 0.1){
                std::cout << sigmanewxy(i) << " " << 0 << std::endl;
            } */
        }
    }
    
    std::cout << "max error on xy is " << maxeps*100.0 << "%" << std::endl;

#endif

#ifdef OUT_SVTU
    std::ofstream out("data.svtu");
    
    double minx=0,miny=0,maxx=0,maxy=0;
    for(int i = 0; i < numnodes; i++){
        minx = std::min(minx,U(2*i));
        miny = std::min(miny,U(2*i + 1));
        maxx = std::max(maxx,U(2*i));
        maxy = std::max(maxy,U(2*i + 1));
    }
    std::cout << "minx: " << minx << " maxx: " << maxx << std::endl;
    std::cout << "miny: " << miny << " maxy: " << maxy << std::endl;

    double minsigmaxx=0,maxsigmaxx=0,minsigmayy=0,maxsigmayy=0,minsigmaxy=0,maxsigmaxy=0;
    for(int i = 0; i < numnodes; i++){
        minsigmaxx = std::min(minsigmaxx,sigmanewxx(i));
        maxsigmaxx = std::max(maxsigmaxx,sigmanewxx(i));

        minsigmayy = std::min(minsigmayy,sigmanewyy(i));
        maxsigmayy = std::max(maxsigmayy,sigmanewyy(i));

        minsigmaxy = std::min(minsigmaxy,sigmanewxy(i));
        maxsigmaxy = std::max(maxsigmaxy,sigmanewxy(i));
    }
    std::cout << "minsigmaxx: " << minsigmaxx << " maxsigmaxx: " << maxsigmaxx << std::endl;
    std::cout << "minsigmayy: " << minsigmayy << " maxsigmayy: " << maxsigmayy << std::endl;
    std::cout << "minsigmaxy: " << minsigmaxy << " maxsigmaxy: " << maxsigmaxy << std::endl;
    //b = color <= 127 ? 255 - color * 2 : 0;
    //g = color <= 127 ? color * 2 : 255 - (color - 127) * 2
    //r = color <= 127 ? 0 : (color - 127) * 2;
    for(int i = 0; i < numnodes; i++){
        double colorx = std::abs(U(2*i) - minx) / std::abs(maxx - minx);
        double colory = std::abs(U(2*i+1) - miny) / std::abs(maxy - miny);
        double colorsigmaxx = std::abs(sigmanewxx(i) - minsigmaxx) / std::abs(maxsigmaxx - minsigmaxx);
        double colorsigmayy = std::abs(sigmanewyy(i) - minsigmayy) / std::abs(maxsigmayy - minsigmayy);
        double colorsigmaxy = std::abs(sigmanewxy(i) - minsigmaxy) / std::abs(maxsigmaxy - minsigmaxy);
        out << colorx << " " << colory << " " << colorsigmaxx << " " << colorsigmayy << " " << colorsigmaxy << std::endl;
    }

    out.close();

#ifdef OUT_COMPARE

    std::vector<std::pair<double,double>> outputDataxx;
    std::vector<std::pair<double,double>> outputDatayy;
    std::vector<std::pair<double,double>> outputDataxy;

    for(int i = 0; i < numnodes; i++){
        if(nodes[i](1) == 0){
            outputDataxx.push_back(std::make_pair(nodes[i](0),sigmanewxx(i)));
            outputDatayy.push_back(std::make_pair(nodes[i](0),sigmanewyy(i)));
            outputDataxy.push_back(std::make_pair(nodes[i](0),sigmanewxy(i)));
            //out << nodes[i](0) << " " << sigmanewxx(i) <<  std::endl;//" " << analyt << std::endl;
        }
    }
    
    std::sort(outputDataxx.begin(),outputDataxx.end());
    std::sort(outputDatayy.begin(),outputDatayy.end());
    std::sort(outputDataxy.begin(),outputDataxy.end());
    out.open("sigmaxx.txt");
    for(auto& i : outputDataxx){
        out << i.first << " " << i.second << std::endl;
    }
    out.close();
    out.open("sigmayy.txt");
    for(auto& i : outputDatayy){
        out << i.first << " " << i.second << std::endl;
    }
    out.close();
    out.open("sigmaxy.txt");
    for(auto& i : outputDataxy){
        out << i.first << " " << i.second << std::endl;
    }

    out.close();

    out.open("analytxx.txt");

    for(auto& i : outputDataxx){
        double analyt = 1.5 * p * (std::pow(1.0 / (i.first),2) - std::pow(1.0 / (i.first),4));
        out << i.first << " " << analyt << std::endl;
    }

     out.close();

     out.open("analytyy.txt");

    for(auto& i : outputDataxx){
        double analyt = 0.5 * p * (2 + std::pow(1.0 / (i.first),2) + 3 * std::pow(1.0 / (i.first),4));
        out << i.first << " " << analyt << std::endl;
    }

     out.close();

     out.open("analytxy.txt");

    for(auto& i : outputDataxx){
        double analyt = 0;
        out << i.first << " " << analyt << std::endl;
    }

     out.close();
#endif
#ifdef SMART_COMPARE

    std::ofstream outxx("sigmaxx.txt");
    std::ofstream axx("analytxx.txt");

    Eigen::Vector2d e(49,0);
    e /= 1000;
    counter = 0;
    for(int i = 1; i <= 1000; i++){
        Eigen::Vector2d point = e * i + Eigen::Vector2d(1,0);
        double valxx = 0;

        for(int j = 0; j < numelements; j++){
            /* Eigen::Vector2d a = nodes[elements[j][1]] - nodes[elements[j][0]];
            Eigen::Vector2d b = nodes[elements[j][2]] - nodes[elements[j][0]];

            Eigen::Matrix2d jacob; 
            jacob << a(0), b(0), a(1), b(1);

            Eigen::Vector2d tmpP = jacob.inverse() * (point - nodes[elements[j][0]]);
            if((tmpP(0) - 1) <= 1e-10 && tmpP(0) >= 1e-10 && (tmpP(1) - 1 - tmpP(0)) <= 1e-10 && tmpP(1) >= 1e-10){
                counter++;
                for(int k = 0; k < 3; k++){
                    valxx += N[3 * j + k].dot(Eigen::Vector3d(1,point(0),point(1))) * sigmanewxx(elements[j][k]);
                }
                break;
            } */
            Eigen::Vector2d a = nodes[elements[j][0]];
            Eigen::Vector2d b = nodes[elements[j][1]];
            Eigen::Vector2d c = nodes[elements[j][2]];

            Eigen::Vector3d ab(b(0) - a(0),b(1) - a(1), 0);
            Eigen::Vector3d bc(c(0) - b(0),c(1) - b(1), 0);
            Eigen::Vector3d ca(a(0) - c(0),a(1) - c(1), 0);

            Eigen::Vector3d ap(point(0) - a(0), point(1) - a(1), 0);
            Eigen::Vector3d bp(point(0) - b(0), point(1) - b(1), 0);
            Eigen::Vector3d cp(point(0) - c(0), point(1) - c(1), 0);

            double v1 = (ab.cross(ap))(2);
            double v2 = (bc.cross(bp))(2);
            double v3 = (ca.cross(cp))(2);

            if((v1 <= 0 && v2 <= 0 && v3 <= 0) || (v1 >= 0 && v2 >= 0 && v3 >= 0)){
                for(int k = 0; k < 3; k++){
                    valxx += N[3 * j + k].dot(Eigen::Vector3d(1,point(0),point(1))) * sigmanewxx(elements[j][k]);
                }
            }

        }

        double analyt = 1.5 * p * (std::pow(1.0 / point(0),2) - std::pow(1.0 / point(0),4));

        outxx << point(0) << " " << valxx << std::endl;
        axx << point(0) << " " << analyt << std::endl;
    }

    outxx.close();
    axx.close();
#endif
#ifdef OUT_BIS
    //out.open("graph.txt");
    outxx.open("graphxx.txt");
    std::ofstream outyy("graphyy.txt");
    std::ofstream outxy("graphxy.txt");

    e = {20,20};
    e /= 500;
    counter = 0;
    for(int i = 1; i <= 500; i++){
        Eigen::Vector2d point = e * i;
        double valxx = 0;
        double valyy = 0;
        double valxy = 0;

        for(int j = 0; j < numelements; j++){
            /* Eigen::Vector2d a = nodes[elements[j][1]] - nodes[elements[j][0]];
            Eigen::Vector2d b = nodes[elements[j][2]] - nodes[elements[j][0]];

            Eigen::Matrix2d jacob; 
            jacob << a(0), b(0), a(1), b(1);

            Eigen::Vector2d tmpP = jacob.inverse() * (point - nodes[elements[j][0]]);
            if((tmpP(0) - 1) <= 1e-10 && tmpP(0) >= 1e-10 && (tmpP(1) - 1 - tmpP(0)) <= 1e-10 && tmpP(1) >= 1e-10){
                counter++;
                for(int k = 0; k < 3; k++){
                    valxx += N[3 * j + k].dot(Eigen::Vector3d(1,point(0),point(1))) * sigmanewxx(elements[j][k]);
                    valyy += N[3 * j + k].dot(Eigen::Vector3d(1,point(0),point(1))) * sigmanewyy(elements[j][k]);
                    valxy += N[3 * j + k].dot(Eigen::Vector3d(1,point(0),point(1))) * sigmanewxy(elements[j][k]);
                }
                break;
            } */

            Eigen::Vector2d a = nodes[elements[j][0]];
            Eigen::Vector2d b = nodes[elements[j][1]];
            Eigen::Vector2d c = nodes[elements[j][2]];

            Eigen::Vector3d ab(b(0) - a(0),b(1) - a(1), 0);
            Eigen::Vector3d bc(c(0) - b(0),c(1) - b(1), 0);
            Eigen::Vector3d ca(a(0) - c(0),a(1) - c(1), 0);

            Eigen::Vector3d ap(point(0) - a(0), point(1) - a(1), 0);
            Eigen::Vector3d bp(point(0) - b(0), point(1) - b(1), 0);
            Eigen::Vector3d cp(point(0) - c(0), point(1) - c(1), 0);

            double v1 = (ab.cross(ap))(2);
            double v2 = (bc.cross(bp))(2);
            double v3 = (ca.cross(cp))(2);

            if((v1 <= 0 && v2 <= 0 && v3 <= 0) || (v1 >= 0 && v2 >= 0 && v3 >= 0)){
                for(int k = 0; k < 3; k++){
                    valxx += N[3 * j + k].dot(Eigen::Vector3d(1,point(0),point(1))) * sigmanewxx(elements[j][k]);
                    valyy += N[3 * j + k].dot(Eigen::Vector3d(1,point(0),point(1))) * sigmanewyy(elements[j][k]);
                    valxy += N[3 * j + k].dot(Eigen::Vector3d(1,point(0),point(1))) * sigmanewxy(elements[j][k]);

                }
            }
        }

        outxx << point.norm() << " " << valxx << std::endl;
        outyy << point.norm() << " " << valyy << std::endl;
        outxy << point.norm() << " " << valxy << std::endl;
    }
    std::cout << counter << std::endl;
    outxx.close();
    outyy.close();
    outxy.close();
    //out.close();
#endif
#endif

    return 0;
}