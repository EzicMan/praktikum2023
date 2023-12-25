#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

//constexpr int numnodes    = 13900;
//constexpr int numelements = 27307;

//constexpr int numnodes    = 2366;
//constexpr int numelements = 4185;

constexpr double p  = -0.01;
constexpr double E  = 1;
constexpr double nu = 0.25;

int main(int argc, char** argv){
    if(argc < 3){
        std::cout << "Not enough parameteres. Usage: " << argv[0] << " <mesh file> <numpoints on graph>" << std::endl;
        return 0;
    }
    std::ifstream in;
    std::vector<std::size_t> fixNodesX;
    std::vector<std::size_t> fixNodesY;
    in.open(argv[1]);
    std::cout << "Opened file. Starting string count" << std::endl;

    for(int i = 0; i < 2; i++){
        in.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    std::string reading = "";
    std::getline(in,reading);
    int numnodes = 0;
    while(reading[0] != '$'){
        numnodes++;
        std::getline(in,reading);
    }

    for(int i = 0; i < 3; i++){
        in.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }

    std::getline(in,reading);
    int numelements = 0;
    while(reading[0] == ' '){
        numelements++;
        std::getline(in,reading);
    }
    
    int numloads = 0;
    while(reading.find("*BOUNDARY_SPC_SET") == std::string::npos){
        if(reading.find("*LOAD_SEGMENT") != std::string::npos){
            numloads++;
        }
        std::getline(in,reading);
    }

    for(int i = 0; i < 3; i++){
        in.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }

    int temporary;
    while(in >> temporary){
        fixNodesX.push_back(temporary-1);
        fixNodesY.push_back(temporary-1);
    }
    
    in.close();

    std::cout << "numnodes: " << numnodes << " numelements:" << numelements << " numloads: " << numloads << std::endl;


    std::vector<Eigen::Vector2d>          nodes(numnodes);
    std::vector<Eigen::Vector2d>          pOnNodes(numnodes, {0,0});
    std::vector<std::vector<std::size_t>> elements(numelements);

    //in.open("../Data_0.5.k");

    std::cout << "Reading data from file" << std::endl;
    in.open(argv[1]);
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
    // for(int i = 0; i < numnodes; i++){
    //     if     (std::abs(nodes[i](0)) < std::numeric_limits<double>::min()) fixNodesX.push_back(i);
    //     else if(std::abs(nodes[i](1)) < std::numeric_limits<double>::min()) fixNodesY.push_back(i);

    //     if(std::abs(nodes[i](1) - ymax) < std::numeric_limits<double>::min()) {
    //         pOnNodes[i](1) = p;
    //     }
    // }

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
    in.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

    for(int i = 0; i < numloads; i++){
        in.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
        double temp;
        int node1, node2;
        in >> temp >> temp >> temp >> node1 >> node2;
        pOnNodes[node1-1](1) = p;
        pOnNodes[node2-1](1) = p;
        for(int i = 0; i < 4; i++){
            in.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
        }
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
        Floc(0) = p;
        Floc(1) = p;
        Floc(2) = p; 
        Floc(3) = p; 
        Floc(4) = p; 
        Floc(5) = p;
        Eigen::Vector2d p1;
        Eigen::Vector2d p2;
        Eigen::Vector2d p3;
        if(n1 && n2){
            p1 = nodes[elements[i][0]];
            p2 = nodes[elements[i][1]];
            p3 = nodes[elements[i][2]];
        }else if(n1 && n3){
            p1 = nodes[elements[i][0]];
            p2 = nodes[elements[i][2]];
            p3 = nodes[elements[i][1]];
        }else if(n2 && n3){
            p1 = nodes[elements[i][1]];
            p2 = nodes[elements[i][2]];
            p3 = nodes[elements[i][0]];
        }

        //for(int j = 0; j < 3; j++){
        //    Floc(2*j+1) *= (ymax * N[3 * i + j](2) + N[3 * i + j](0)) * std::abs(xi1 - xi) +
        //    N[3 * i + j](1) * 0.5 * std::abs(xi1 - xi) * (xi1 + xi);
        //}
        if(p1(0) > p2(0)){
            std::swap(p1,p2);
        }
        Eigen::Vector2d normal({(p1-p2)(1),-(p1-p2)(0)});
        if(normal.dot(p3-p1) > 0){
            normal = -normal;
        }

        normal.normalize();
        
        double k2 = (p2(1) - p1(1))/(p2(0) - p1(0));
        double k1 = p1(1) - k2 * p1(0);

        for(int j = 0; j < 3; j++){
            double a = N[3*i+j](0);
            double b = N[3*i+j](1);
            double c = N[3*i+j](2);
            Floc(2*j)   *= -0.5 * std::sqrt(k2*k2 + 1) * (p1(0) - p2(0)) * (2*a+b*(p1(0) + p2(0)) + 2 * c * k1 + c * k2 * (p1(0) + p2(0))) * normal(0);
            Floc(2*j+1) *= -0.5 * std::sqrt(k2*k2 + 1) * (p1(0) - p2(0)) * (2*a+b*(p1(0) + p2(0)) + 2 * c * k1 + c * k2 * (p1(0) + p2(0))) * normal(1);
        }

        
        Fglob.coeffRef(2*elements[i][0]    ) += Floc(0);
        Fglob.coeffRef(2*elements[i][0] + 1) += Floc(1);
        Fglob.coeffRef(2*elements[i][1]    ) += Floc(2);
        Fglob.coeffRef(2*elements[i][1] + 1) += Floc(3);
        Fglob.coeffRef(2*elements[i][2]    ) += Floc(4);
        Fglob.coeffRef(2*elements[i][2] + 1) += Floc(5);
    }

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
    /* double minsigmaxx=0,maxsigmaxx=0,minsigmayy=0,maxsigmayy=0,minsigmaxy=0,maxsigmaxy=0;
    for(int i = 0; i < numnodes; i++){
        minsigmaxx = std::min(minsigmaxx,sigma[i](0));
        maxsigmaxx = std::max(maxsigmaxx,sigma[i](0));

        minsigmayy = std::min(minsigmayy,sigma[i](1));
        maxsigmayy = std::max(maxsigmayy,sigma[i](1));

        minsigmaxy = std::min(minsigmaxy,sigma[i](2));
        maxsigmaxy = std::max(maxsigmaxy,sigma[i](2));
    }
    std::cout << "minsigmaxx: " << minsigmaxx << " maxsigmaxx: " << maxsigmaxx << std::endl;
    std::cout << "minsigmayy: " << minsigmayy << " maxsigmayy: " << maxsigmayy << std::endl;
    std::cout << "minsigmaxy: " << minsigmaxy << " maxsigmaxy: " << maxsigmaxy << std::endl; */


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
                double b1tmp = b1;
                b1 = b1 * jacob(0,0) + c1 * jacob(1,0);
                c1 = b1tmp * jacob(0,1) + c1 * jacob(1,1);

                double a2 = N[3 * i + k](0);
                double b2 = N[3 * i + k](1);
                double c2 = N[3 * i + k](2);

                a2 = a2 + b2 * nodes[elements[i][0]](0) + c2 * nodes[elements[i][0]](1);
                double b2tmp = b2;
                b2 = b2 * jacob(0,0) + c2 * jacob(1,0);
                c2 = b2tmp * jacob(0,1) + c2 * jacob(1,1);

                Cloc(j,k) = jacob.determinant() / 24.0 * (4 * a1 * (3 * a2 + b2 + c2) + 4 * a2 * (b1 + c1) + 2 * b1 * b2 + b1 * c2 + b2 * c1 + 2 * c1 * c2);
            }
        }
        for(int j = 0; j < 3; j++){
            double a1 = N[3 * i + j](0);
            double b1 = N[3 * i + j](1);
            double c1 = N[3 * i + j](2);

            a1 = a1 + b1 * nodes[elements[i][0]](0) + c1 * nodes[elements[i][0]](1);
            double b1tmp = b1;
            b1 = b1 * jacob(0,0) + c1 * jacob(1,0);
            c1 = b1tmp * jacob(0,1) + c1 * jacob(1,1);

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

    Eigen::VectorXd sigmanewxx = ldlt.solve(Rglobxx);
    Eigen::VectorXd sigmanewyy = ldlt.solve(Rglobyy);
    Eigen::VectorXd sigmanewxy = ldlt.solve(Rglobxy);

    Eigen::Vector3d sigmaLocal({0,0,0});
    constexpr double fullvolume = 3.1415926 * 100;
    double volume = 0;
    for(int i = 0; i < numelements; i++){
        Eigen::Vector2d a = nodes[elements[i][1]] - nodes[elements[i][0]];
        Eigen::Vector2d b = nodes[elements[i][2]] - nodes[elements[i][0]];

        Eigen::Matrix2d jacob; 
        jacob << a(0), b(0), a(1), b(1);

        volume += jacob.determinant() * 0.5;

        sigmaLocal(0) += sigma[i](0);
        sigmaLocal(1) += sigma[i](1);
        sigmaLocal(2) += sigma[i](2);
    }
    sigmaLocal /= fullvolume;
    double porosity = (fullvolume - volume) / fullvolume;
    std::cout << (sigmaLocal(0) / p + porosity) << " " << (sigmaLocal(1) / p + porosity) << " " << (sigmaLocal(2) / p) << std::endl;
    std::cout << porosity << std::endl;
    /* double maxeps = 0;

    for(int i = 0; i < numnodes; i++){
        if(nodes[i](1) == 0){
            double analyt = 0.5 * p * (2 + std::pow(1.0 / (nodes[i](0)),2) + 3 * std::pow(1.0 / (nodes[i](0)),4));
            maxeps = std::max(maxeps, std::abs((sigmanewyy(i) - analyt) / analyt));
        }
    }
    
    std::cout << "max error on yy is " << maxeps*100.0 << "%" << std::endl;

    maxeps = 0;

    for(int i = 0; i < numnodes; i++){
        if(nodes[i](1) == 0){
            double analyt = 1.5 * p * (std::pow(1.0 / (nodes[i](0)),2) - std::pow(1.0 / (nodes[i](0)),4));
            if(analyt == 0){
                maxeps = std::max(maxeps, std::abs(sigmanewxx(i)));
            }else{
                maxeps = std::max(maxeps, std::abs((sigmanewxx(i) - analyt) / analyt));
            }
        }
    }
    
    std::cout << "max error on xx is " << maxeps*100.0 << "%" << std::endl;

    maxeps = 0;

    for(int i = 0; i < numnodes; i++){
        if(nodes[i](1) == 0){
            maxeps = std::max(maxeps, std::abs(sigmanewxy(i)));
        }
    }
    
    std::cout << "max error on xy is " << maxeps*100.0 << "%" << std::endl; */

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
    for(int i = 0; i < numnodes; i++){
        double colorx = std::abs(U(2*i) - minx) / std::abs(maxx - minx);
        double colory = std::abs(U(2*i+1) - miny) / std::abs(maxy - miny);
        double colorsigmaxx, colorsigmayy, colorsigmaxy;
        colorsigmaxx = std::abs(sigmanewxx(i) - minsigmaxx) / std::abs(maxsigmaxx - minsigmaxx);
        colorsigmayy = std::abs(sigmanewyy(i) - minsigmayy) / std::abs(maxsigmayy - minsigmayy);
        colorsigmaxy = std::abs(sigmanewxy(i) - minsigmaxy) / std::abs(maxsigmaxy - minsigmaxy);
        out << colorx << " " << colory << " " << colorsigmaxx << " " << colorsigmayy << " " << colorsigmaxy << std::endl;
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
    out.open(std::string(argv[1])+"_graph.txt");
    //std::ofstream outyy("graphyy.txt");
    //std::ofstream outxy("graphxy.txt");

    Eigen::Vector2d s = {-10,-10};
    Eigen::Vector2d e = {10,10};
    int numpoints = atoi(argv[2]);
    Eigen::Vector2d sm = e - s;
    sm /= numpoints;
    for(int i = 1; i <= numpoints; i++){
        Eigen::Vector2d point = sm * i + s;
        double valxx = 0;
        double valyy = 0;
        double valxy = 0;

        for(int j = 0; j < numelements; j++){
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

        out << (point-s).norm() << " " << valxx << " " << valyy << " " << valxy << std::endl;
        //outyy << point.norm() << " " << valyy << std::endl;
        //outxy << point.norm() << " " << valxy << std::endl;
    }
    out.close();
    //outyy.close();
    //outxy.close();
#endif

    return 0;
}