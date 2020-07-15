float calculateLocalD(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, YE, 2, 1, m));
    row1.push_back(calcularTenedor(e, ZETA, 2, 1, m));

    row2.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, ZETA, 3, 1, m));

    row3.push_back(calcularTenedor(e, EQUIS, 4, 1, m));
    row3.push_back(calcularTenedor(e, YE, 4, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

void calculateAlpha(int i,Matrix &A,mesh m){
    zeroes(A,3);
    element e = m.getElement(i);
    
    A.at(0).at(0) = OperarRestaTenedor(e, YE, ZETA, 3, 4, m);
    A.at(0).at(1) = OperarRestaTenedor(e, YE, ZETA, 4, 2, m);
    A.at(0).at(2) = OperarRestaTenedor(e, YE, ZETA, 2, 3, m);

    A.at(1).at(0) = OperarRestaTenedor(e, EQUIS, ZETA, 4, 3, m);
    A.at(1).at(1) = OperarRestaTenedor(e, EQUIS, ZETA, 2, 4, m);
    A.at(1).at(2) = OperarRestaTenedor(e, EQUIS, ZETA, 3, 2, m);

    A.at(2).at(0) = OperarRestaTenedor(e, EQUIS, YE, 3, 4, m);
    A.at(2).at(1) = OperarRestaTenedor(e, EQUIS, YE, 4, 2, m);
    A.at(2).at(2) = OperarRestaTenedor(e, EQUIS, YE, 2, 3, m);

}

void calculateBeta(Matrix &B){
    zeroes(B,3,12);

    B.at(0).at(0) = -1; 
    B.at(0).at(1) =  1; 
    B.at(0).at(2) =  0; 
    B.at(0).at(3) =  0; 
    B.at(0).at(4) = -1; 
    B.at(0).at(5) =  1; 
    B.at(0).at(6) =  0; 
    B.at(0).at(7) =  0; 
    B.at(0).at(8) = -1; 
    B.at(0).at(9) =  1; 
    B.at(0).at(10) =  0; 
    B.at(0).at(11) =  0; 
    
    B.at(1).at(0) = -1; 
    B.at(1).at(1) = 0; 
    B.at(1).at(2) = 1;
    B.at(1).at(3) = 0;
    B.at(1).at(4) = -1; 
    B.at(1).at(5) = 0; 
    B.at(1).at(6) = 1;
    B.at(1).at(7) = 0;
    B.at(1).at(8) = -1; 
    B.at(1).at(9) = 0; 
    B.at(1).at(10) = 1;
    B.at(1).at(11) = 0;

    B.at(2).at(0) = -1;
    B.at(2).at(1) = 0;
    B.at(2).at(2) = 0;
    B.at(2).at(3) = 1;
    B.at(2).at(4) = -1;
    B.at(2).at(5) = 0;
    B.at(2).at(6) = 0;
    B.at(2).at(7) = 1;
    B.at(2).at(8) = -1;
    B.at(2).at(9) = 0;
    B.at(2).at(10) = 0;
    B.at(2).at(11) = 1;
}

float calculateBK(mesh m, int i){
	element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float z1, z2, z3, z4;
    z1 = n1.getZ();
    z2 = n2.getZ();
    z3 = n3.getZ();
    z4 = n4.getZ();

    float y1, y2, y3, y4;
    y1 = n1.getY();
    y2 = n2.getY();
    y3 = n3.getY();
    y4 = n4.getY();

    float x1, x2, x3, x4;
    x1 = n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();
    
    if(z1<0){
    	z1=z1*(-1);
	};
	  if(z2<0){
    	z2=z2*(-1);
	};
	  if(z3<0){
    	z3=z3*(-1);
	};
	  if(z4<0){
    	z4=z4*(-1);
	};
	
	  if(y1<0){
    	y1=y1*(-1);
	};
	  if(y2<0){
    	y2=y2*(-1);
	};
	  if(y3<0){
    	y3=y3*(-1);
	};
	  if(y4<0){
    	y4=y4*(-1);
	};
	
	  if(x1<0){
    	x1=x1*(-1);
	};
	  if(x2<0){
    	x2=x2*(-1);
	};
	  if(x3<0){
    	x3=x3*(-1);
	};
	  if(z4<0){
    	x4=x4*(-1);
	};

	if(z1==0){
    	z1=1;
	};
	  if(z2==0){
    	z2=1;
	};
	  if(z3==0){
    	z3=1;
	};
	  if(z4==0){
    	z4=1;
	};
	
	  if(y1==0){
    	y1=1;
	};
	  if(y2==0){
    	y2=1;
	};
	  if(y3==0){
    	y3=1;
	};
	  if(y4==0){
    	y4=1;
	};
	
	  if(x1==0){
    	x1=1;
	};
	  if(x2==0){
    	x2=1;
	};
	  if(x3==0){
    	x3=1;
	};
	  if(x4==0){
    	x4=1;
	};
 

    

    //return ((8*(sqrt(x2+10)-sqrt(x3+10)) *pow(x3,2)/(105*(x2-x3)*(x3-x4)))-((8*(sqrt(x2+10)-sqrt(x4+10))* pow(x4,2))/(105*(x2-x4)*(x3-x4)))+((64* pow(x2,2.5)+35*x3*x4*(y2+y3+y4+z2+z3+z4))/(840*x3*x4))+(( pow(x2,2.5)*(((8)/(105*x3))-((8)/(105*x4))))/(x3-x4))+((sqrt(x2+10)*(((8*x3)/(105))-((8*x4)/(105))))/(x3-x4)));
    return 12;
    
}



void calculateOmega(Matrix &C){
    zeroes(C,3,4);
    C.at(0).at(0) = -1; C.at(0).at(1) = 1; C.at(0).at(2) = 0; C.at(0).at(3) = 0;
    C.at(1).at(0) = -1; C.at(1).at(1) = 0; C.at(1).at(2) = 1; C.at(1).at(3) = 0;
    C.at(2).at(0) = -1; C.at(2).at(1) = 0; C.at(2).at(2) = 0; C.at(2).at(3) = 1;
}

void calculateGamma(Matrix &m){
	zeroes(m,12,3);

	m.at(0).at(0) = 1;   m.at(0).at(1) = 0;   m.at(0).at(2) = 0;
	m.at(1).at(0) = 1;   m.at(1).at(1) = 0;   m.at(1).at(2) = 0; 
    m.at(2).at(0) = 1;   m.at(2).at(1) = 0;   m.at(2).at(2) = 0;
	m.at(3).at(0) = 1;   m.at(3).at(1) = 0;   m.at(3).at(2) = 0; 
    m.at(4).at(0) = 0;   m.at(4).at(1) = 1;   m.at(4).at(2) = 0;
	m.at(5).at(0) = 0;   m.at(5).at(1) = 1;   m.at(5).at(2) = 0; 
    m.at(6).at(0) = 0;   m.at(6).at(1) = 1;   m.at(6).at(2) = 0;
	m.at(7).at(0) = 0;   m.at(7).at(1) = 1;   m.at(7).at(2) = 0; 
    m.at(8).at(0) = 0;   m.at(8).at(1) = 0;   m.at(8).at(2) = 1;
	m.at(9).at(0) = 0;   m.at(9).at(1) = 0;   m.at(9).at(2) = 1; 
    m.at(10).at(0) = 0;  m.at(10).at(1) = 0;  m.at(10).at(2) = 1;
	m.at(11).at(0) = 0;  m.at(11).at(1) = 0;  m.at(11).at(2) = 1; 
	
}
//MATRIZ A
void calculateGammaA(Matrix &G1, mesh m, int i){
	zeroes(G1,12,3);
	
	element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1());
    node n2 = m.getNode(e.getNode2());
    node n3 = m.getNode(e.getNode3());
    node n4 = m.getNode(e.getNode4());

    float z1, z2, z3, z4;
    z1 = n1.getZ();
    z2 = n2.getZ();
    z3 = n3.getZ();
    z4 = n4.getZ();

    float y1, y2, y3, y4;
    y1 = n1.getY();
    y2 = n2.getY();
    y3 = n3.getY();
    y4 = n4.getY();

    float x1, x2, x3, x4;
    x1 = n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();

    if(z1<0){
    	z1=z1*(-1);
	};
	  if(z2<0){
    	z2=z2*(-1);
	};
	  if(z3<0){
    	z3=z3*(-1);
	};
	  if(z4<0){
    	z4=z4*(-1);
	};
	
	  if(y1<0){
    	y1=y1*(-1);
	};
	  if(y2<0){
    	y2=y2*(-1);
	};
	  if(y3<0){
    	y3=y3*(-1);
	};
	  if(y4<0){
    	y4=y4*(-1);
	};
	
	  if(x1<0){
    	x1=x1*(-1);
	};
	  if(x2<0){
    	x2=x2*(-1);
	};
	  if(x3<0){
    	x3=x3*(-1);
	};
	  if(z4<0){
    	x4=x4*(-1);
	};


 	if(z1==0){
    	z1=1;
	};
	  if(z2==0){
    	z2=1;
	};
	  if(z3==0){
    	z3=1;
	};
	  if(z4==0){
    	z4=1;
	};
	
	  if(y1==0){
    	y1=1;
	};
	  if(y2==0){
    	y2=1;
	};
	  if(y3==0){
    	y3=1;
	};
	  if(y4==0){
    	y4=1;
	};
	
	  if(x1==0){
    	x1=1;
	};
	  if(x2==0){
    	x2=1;
	};
	  if(x3==0){
    	x3=1;
	};
	  if(x4==0){
    	x4=1;
	};



    float a = ((16*(sqrt(x2+10)-sqrt(x3+10))*pow(x3,2))/(945*(x2-x3)*(x3-x4)))-((16*(sqrt(x2+10)-sqrt(x4+10))* pow(x4,2))/(945*(x2-x4)*(x3-x4)))+((64* pow(x2,3.5))*(x3+x4)+64* pow(x2,2.5)*x3*x4+63*pow(x3,2)*pow(x4,2)*(y2+y3+y4))/(3780*pow(x3,2)*pow(x4,2))+(( pow(x2,3.5)*(((16)/(945*pow(x3,2)))-((16)/(945*pow(x4,2)))))/(x3-x4))+(( pow(x2,2.5) *(((16)/(945*x3))-((16)/(945*x4))))/(x3-x4))+((sqrt(x2+10)*(((16*x3)/(945))-((16*x4)/(945))))/(x3-x4));
    float b = ((32*pow(x2,2.5)+9*x3*x4*(2*y2+y3+y4))/(540*x3*x4))+((8*( pow(x2,1.5)-3*sqrt(x2+10)*x3+2* pow(x3,1.5)*pow(x2,2)/(945*pow(x2-x3,2)*(x3-x4)))-((8*(pow(x2,1.5)-3*sqrt(x2+10)*x4+2* pow(x4,1.5))*pow(x4,2))/(945*pow(x2-x4,2)*(x3-x4)))+( pow(x2,2.5) *(((8)/(135*x3))-((8)/(135*x4))))/(x3-x4))+((sqrt(x2+10)*(((8*x3)/(315))-((8*x4)/(315))))/(x3-x4)));
    float c = ((16*(sqrt(x2+10)-sqrt(x3+10))*pow(x3,2))/(945*(x2-x3)*(x3-x4)))-((16*(sqrt(x2+10)-sqrt(x4+10))* pow(x4,2))/(945*(x2-x4)*(x3-x4)))+((64* pow(x2,3.5))*(x3+x4)+64* pow(x2,2.5)*x3*x4+63*pow(x3,2)*pow(x4,2)*(y2+y3+y4))/(3780*pow(x3,2)*pow(x4,2))+(( pow(x2,3.5)*(((16)/(945*pow(x3,2)))-((16)/(945*pow(x4,2)))))/(x3-x4))+(( pow(x2,2.5) *(((16)/(945*x3))-((16)/(945*x4))))/(x3-x4))+((sqrt(x2+10)*(((16*x3)/(945))-((16*x4)/(945))))/(x3-x4));
    float d = ((32*pow(x2,2.5)+9*x3*x4*(2*y2+y3+y4))/(540*x3*x4))+((8*( pow(x2,1.5)-3*sqrt(x2+10)*x3+2* pow(x3,1.5)*pow(x2,2)/(945*pow(x2-x3,2)*(x3-x4)))-((8*(pow(x2,1.5)-3*sqrt(x2+10)*x4+2* pow(x4,1.5))*pow(x4,2))/(945*pow(x2-x4,2)*(x3-x4)))+( pow(x2,2.5) *(((8)/(135*x3))-((8)/(135*x4))))/(x3-x4))+((sqrt(x2+10)*(((8*x3)/(315))-((8*x4)/(315))))/(x3-x4)));
 
 	a=2;
 	b=3;
 	c=4;
 	d=5;

	G1.at(0).at(0) = a;                                     G1.at(0).at(1) = 0;   									G1.at(0).at(2) = 0;
	G1.at(1).at(0) = b;                                     G1.at(1).at(1) = 0;   									G1.at(1).at(2) = 0; 
    G1.at(2).at(0) = c;                                     G1.at(2).at(1) = 0;   									G1.at(2).at(2) = 0;
	G1.at(3).at(0) = d;                                     G1.at(3).at(1) = 0;   									G1.at(3).at(2) = 0; 
    G1.at(4).at(0) = 0;   									G1.at(4).at(1) = a;                                     G1.at(4).at(2) = 0;
	G1.at(5).at(0) = 0;   									G1.at(5).at(1) = b;                                     G1.at(5).at(2) = 0; 
    G1.at(6).at(0) = 0;   									G1.at(6).at(1) = c;                                     G1.at(6).at(2) = 0;
	G1.at(7).at(0) = 0;   									G1.at(7).at(1) = d;                                     G1.at(7).at(2) = 0; 
    G1.at(8).at(0) = 0;   									G1.at(8).at(1) = 0;   									G1.at(8).at(2) = a; 
	G1.at(9).at(0) = 0;   									G1.at(9).at(1) = 0;   									G1.at(9).at(2) = b; 
    G1.at(10).at(0) = 0; 									G1.at(10).at(1) = 0; 									G1.at(10).at(2) = c;
	G1.at(11).at(0) = 0;  									G1.at(11).at(1) = 0;  									G1.at(11).at(2) = d; 
}

//matriz G
void calculateGammaG(Matrix &G1, mesh m, int i){
	zeroes(G1,12,3);
	
	element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float z1, z2, z3, z4;
    z1 = n1.getZ();
    z2 = n2.getZ();
    z3 = n3.getZ();
    z4 = n4.getZ();

    float y1, y2, y3, y4;
    y1 = n1.getY();
    y2 = n2.getY();
    y3 = n3.getY();
    y4 = n4.getY();

    float x1, x2, x3, x4;
    x1 = n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();


 if(z1<0){
    	z1=z1*(-1);
	};
	  if(z2<0){
    	z2=z2*(-1);
	};
	  if(z3<0){
    	z3=z3*(-1);
	};
	  if(z4<0){
    	z4=z4*(-1);
	};
	
	  if(y1<0){
    	y1=y1*(-1);
	};
	  if(y2<0){
    	y2=y2*(-1);
	};
	  if(y3<0){
    	y3=y3*(-1);
	};
	  if(y4<0){
    	y4=y4*(-1);
	};
	
	  if(x1<0){
    	x1=x1*(-1);
	};
	  if(x2<0){
    	x2=x2*(-1);
	};
	  if(x3<0){
    	x3=x3*(-1);
	};
	  if(z4<0){
    	x4=x4*(-1);
	};

if(z1==0){
    	z1=1;
	};
	  if(z2==0){
    	z2=1;
	};
	  if(z3==0){
    	z3=1;
	};
	  if(z4==0){
    	z4=1;
	};
	
	  if(y1==0){
    	y1=1;
	};
	  if(y2==0){
    	y2=1;
	};
	  if(y3==0){
    	y3=1;
	};
	  if(y4==0){
    	y4=1;
	};
	
	  if(x1==0){
    	x1=1;
	};
	  if(x2==0){
    	x2=1;
	};
	  if(x3==0){
    	x3=1;
	};
	  if(x4==0){
    	x4=1;
	};

 	z1=3;
  	z2=3;
   	z3=3;
    z4=3;
    x1=3;
  	x2=3;
   	x3=3;
    x4=3;
    y1=3;
  	y2=3;
  	y3=3;
    y4=3;

    float a = ((y2+y3+y4)/(120));
    float b = ((2*y2+y3+y4)/(120));
    float c = ((y2+2*y3+y4)/(120));
    float d = ((y2+y3+2*y4)/(120));
 

	G1.at(0).at(0) = a;                                     G1.at(0).at(1) = 0;   									G1.at(0).at(2) = 0;
	G1.at(1).at(0) = b;                                     G1.at(1).at(1) = 0;   									G1.at(1).at(2) = 0; 
    G1.at(2).at(0) = c;                                     G1.at(2).at(1) = 0;   									G1.at(2).at(2) = 0;
	G1.at(3).at(0) = d;                                     G1.at(3).at(1) = 0;   									G1.at(3).at(2) = 0; 
    G1.at(4).at(0) = 0;   									G1.at(4).at(1) = a;                                     G1.at(4).at(2) = 0;
	G1.at(5).at(0) = 0;   									G1.at(5).at(1) = b;                                     G1.at(5).at(2) = 0; 
    G1.at(6).at(0) = 0;   									G1.at(6).at(1) = c;                                     G1.at(6).at(2) = 0;
	G1.at(7).at(0) = 0;   									G1.at(7).at(1) = d;                                     G1.at(7).at(2) = 0; 
    G1.at(8).at(0) = 0;   									G1.at(8).at(1) = 0;   									G1.at(8).at(2) = a; 
	G1.at(9).at(0) = 0;   									G1.at(9).at(1) = 0;   									G1.at(9).at(2) = b; 
    G1.at(10).at(0) = 0; 									G1.at(10).at(1) = 0; 									G1.at(10).at(2) = c;
	G1.at(11).at(0) = 0;  									G1.at(11).at(1) = 0;  									G1.at(11).at(2) = d; 
}

//MATRIZ D
void calculateGammaD(Matrix &G1, mesh m, int i){
	zeroes(G1,12,3);
	
	element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    float z1, z2, z3, z4;
    z1 = n1.getZ();
    z2 = n2.getZ();
    z3 = n3.getZ();
    z4 = n4.getZ();

    float y1, y2, y3, y4;
    y1 = n1.getY();
    y2 = n2.getY();
    y3 = n3.getY();
    y4 = n4.getY();

    float x1, x2, x3, x4;
    x1 = n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();
    
   if(z1<0){
    	z1=z1*(-1);
	};
	  if(z2<0){
    	z2=z2*(-1);
	};
	  if(z3<0){
    	z3=z3*(-1);
	};
	  if(z4<0){
    	z4=z4*(-1);
	};
	
	  if(y1<0){
    	y1=y1*(-1);
	};
	  if(y2<0){
    	y2=y2*(-1);
	};
	  if(y3<0){
    	y3=y3*(-1);
	};
	  if(y4<0){
    	y4=y4*(-1);
	};
	
	  if(x1<0){
    	x1=x1*(-1);
	};
	  if(x2<0){
    	x2=x2*(-1);
	};
	  if(x3<0){
    	x3=x3*(-1);
	};
	  if(x4<0){
    	x4=x4*(-1);
	};

 	if(z1==0){
    	z1=1;
	};
	  if(z2==0){
    	z2=1;
	};
	  if(z3==0){
    	z3=1;
	};
	  if(z4==0){
    	z4=1;
	};
	
	  if(y1==0){
    	y1=1;
	};
	  if(y2==0){
    	y2=1;
	};
	  if(y3==0){
    	y3=1;
	};
	  if(y4==0){
    	y4=1;
	};
	
	  if(x1==0){
    	x1=1;
	};
	  if(x2==0){
    	x2=1;
	};
	  if(x3==0){
    	x3=1;
	};
	  if(x4==0){
    	x4=1;
	};
	
	z1=3;
  	z2=3;
   	z3=3;
    z4=3;
    x1=3;
  	x2=3;
   	x3=3;
    x4=3;
    y1=3;
 	y2=3;
   	y3=3;
    y4=3;
	
    float a = (pow(x3,2)/(24*((x2-x3)*(x3-x4))))-(pow(x4,2)/((24*((x2-x4)*(x3-x4))+(0.04)))*log10(x2+10)*((48*y2+48*y3+48*y4-125*log10(x4+10))/(1440)));
    float b = (pow(x4,3)/(24*(pow(x2-x4,2)*(x3-x4))))-( pow(x3,3)/(24*pow(x2-x3,2)-(x3-x4)))*log10(x2+10)-x4*(60*log10(x4+10)+96*y2+48*y3+48*y4-125/(1440));
    float c = (pow(x4,3)/(24*(pow(x2-x4,2)*(x3-x4))+(48*y2+96*y3+48*y4-125*log10(x4+10))/(1440)));
    float d = (pow(x3,3)/(24*(x2-x3))* pow(x2-x3,2)-(((x2*(3*x3-2*x4)-(2*x3-x4)*x4)* pow(x4,2))/(24*(pow(x2-x4,2) *pow(x2-x3,2)) + ((48*y2+48*y3+96*y4-125*log10(x4+10))/(1440)))));
 
 	a=2;
 	b=3;
 	c=4;
 	d=5;

	G1.at(0).at(0) = a;                                     G1.at(0).at(1) = 0;   									G1.at(0).at(2) = 0;
	G1.at(1).at(0) = b;                                     G1.at(1).at(1) = 0;   									G1.at(1).at(2) = 0; 
    G1.at(2).at(0) = c;                                     G1.at(2).at(1) = 0;   									G1.at(2).at(2) = 0;
	G1.at(3).at(0) = d;                                     G1.at(3).at(1) = 0;   									G1.at(3).at(2) = 0; 
    G1.at(4).at(0) = 0;   									G1.at(4).at(1) = a;                                     G1.at(4).at(2) = 0;
	G1.at(5).at(0) = 0;   									G1.at(5).at(1) = b;                                     G1.at(5).at(2) = 0; 
    G1.at(6).at(0) = 0;   									G1.at(6).at(1) = c;                                     G1.at(6).at(2) = 0;
	G1.at(7).at(0) = 0;   									G1.at(7).at(1) = d;                                     G1.at(7).at(2) = 0; 
    G1.at(8).at(0) = 0;   									G1.at(8).at(1) = 0;   									G1.at(8).at(2) = a; 
	G1.at(9).at(0) = 0;   									G1.at(9).at(1) = 0;   									G1.at(9).at(2) = b; 
    G1.at(10).at(0) = 0; 									G1.at(10).at(1) = 0; 									G1.at(10).at(2) = c;
	G1.at(11).at(0) = 0;  									G1.at(11).at(1) = 0;  									G1.at(11).at(2) = d; 
}


float calculateLocalJ(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 4, 1, m));

    row2.push_back(calcularTenedor(e, YE, 2, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 4, 1, m));

    row3.push_back(calcularTenedor(e, ZETA, 2, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 3, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

Matrix createLocalM(int e,mesh &m){
    Matrix matrixA,matrixK,matrixG,matrixD,matrixC,matrixE,matrixW;
    float u_bar,nu,rho,Ve,J,Determinant, Vol;
    
    /* [ A+K  G ]
       [  D   0 ]
    */
    

    //Matrix A
    Matrix g_matrix, Alpha, Beta, gA_matrix,gG_matrix,gD_matrix;

    //u_bar = m.getParameter(ADJECTIVE_VELOCITY);
    Determinant = calculateLocalD(e,m);
    J = calculateLocalJ(e,m);

    if(Determinant == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    float real_a = (float) (J)/(Determinant);
    
    calculateGammaA(gA_matrix,m,e);
    calculateGammaG(gG_matrix,m,e);
    calculateGammaD(gD_matrix,m,e);

    calculateGamma(g_matrix);
    calculateAlpha(e,Alpha,m);
    calculateBeta(Beta);

    productRealMatrix(real_a, productMatrixMatrix(gA_matrix,productMatrixMatrix(Alpha,Beta,3,3,12),12,3,12),matrixA);
    
    
    //Matrix K
    Matrix Alpha_t,Beta_t,Beta_t_z,Beta_t_z2;

    //nu = m.getParameter(DYNAMIC_VISCOSITY);
    Vol = calculateBK(m,e);
    
    float real_k = (float) (1)/(24*Determinant*Determinant);

    transpose(Alpha,Alpha_t);
    transpose(Beta,Beta_t);

    productRealMatrix(real_k,productMatrixMatrix(Beta_t,productMatrixMatrix(Alpha_t,productMatrixMatrix(Alpha,Beta,3,3,12),3,3,12),12,3,12),matrixK);
    
   

	//Matrix G
	Matrix Omega;
	calculateOmega(Omega);
	float real_g = (float) (J)/(Determinant);
	 
	productRealMatrix(real_g, productMatrixMatrix(gG_matrix,productMatrixMatrix(Alpha,Beta,3,3,12),12,3,12),matrixG);
    

    //Matrix D
    Matrix g_matrix_t,Omega_t;
    float real_d = (float)(J/(Determinant));

    transpose(Omega, Omega_t);
    transpose(gD_matrix,g_matrix_t);
    productRealMatrix(real_d,productMatrixMatrix(Omega_t,productMatrixMatrix(Alpha_t,g_matrix_t,3,3,12),4,3,12),matrixD);


    //Matrix M
    Matrix M;
    zeroes(M,16);
    ubicarSubMatriz(M,0,11,0,11, sumMatrix(matrixA,matrixK,12,12));
    ubicarSubMatriz(M,0,11,12,15,matrixG);
    ubicarSubMatriz(M,12,15,0,11,matrixD);

    return M;
}

//CALCULAR F
void calculateF(Vector &f, mesh &m){
    zeroes(f,3);

    f.at(0) += m.getParameter(EXTERNAL_FORCE_X);
    f.at(1) += m.getParameter(EXTERNAL_FORCE_Y);
    f.at(2) += m.getParameter(EXTERNAL_FORCE_Z);

}

//CALCULAR H

void calculateH(Vector &h, mesh &m){
    zeroes(h,3);

    h.at(0) += m.getParameter(EXTERNAL_FORCE_X);
    h.at(1) += m.getParameter(EXTERNAL_FORCE_Y);
    h.at(2) += m.getParameter(EXTERNAL_FORCE_Z);

}

Vector createLocalb(int e,mesh &m){
    float J;
    Vector b,b_aux,f;
    Matrix g_matrix;

    calculateF(f, m);

    calculateGamma(g_matrix);

    J = calculateLocalJ(e,m);

    if(J == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    zeroes(b_aux,24);
    productMatrixVector(g_matrix,f,b_aux);
    productRealVector(J/24,b_aux,b);
    
    return b;
}
