#ifndef DEF_QUADRIC
#define DEF_QUADRIC
 #include <iostream>

#include <armadillo>



class Quadric
{

public:
    
    arma::mat Q ;


    Quadric(arma::vec normale, arma::vec point)
    {
        Q=arma::mat(4,4),arma::fill::ones;
        double a =normale(0);
        double b =normale(1);
        double c =normale(2);
        double d=-arma::dot(normale,point);

        Q(0,0)=a*a; Q(0,1)=a*b; Q(0,2)=a*c; Q(0,3)=a*d;
        
        Q(1,0)=a*b; Q(1,1)=b*b; Q(1,2)=b*c; Q(1,3)=b*d;   
        
        Q(2,0)=a*c; Q(2,1)=b*c; Q(2,2)=c*c; Q(2,3)=c*d;   
     
        Q(3,0)=a*d; Q(3,1)=b*d; Q(3,2)=c*d; Q(3,3)=d*d;
       
       
    }
void add(Quadric Q2)
{
    Q=Q+Q2.Q;
}
   
double error(arma::vec v)
{


  arma::mat vt=  arma::trans(v);
  arma::mat res=(vt*Q)*v;
    return res(0,0);
}
  

} ;
#endif