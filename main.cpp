#include <iostream>
#include <string>
#include <vector>
#include <math.h> 
#include<pmp/io/io.h>
#include <pmp/SurfaceMesh.h>
#include <pmp/Types.h>
#include <pmp/algorithms/Curvature.h>
#include<pmp/BoundingBox.h>
#include<pmp/utilities.h>
#include <pmp/algorithms/Normals.h>
#include <pmp/MatVec.h>
#include "Quadric.h"
#include<pmp/algorithms/Quadric.h>
using namespace std;







void computeQuadrix(pmp::SurfaceMesh & ppmmesh)
{
  pmp::Normals::compute_face_normals(ppmmesh);
  auto quadrics =ppmmesh.add_vertex_property<pmp::Scalar>("added:quadrics");
  auto normales= ppmmesh.get_face_property<pmp::Normal>("f:normal");
  pmp::VertexProperty<pmp::Point>  points = ppmmesh.get_vertex_property<pmp::Point>("v:point");


 //
  //cout<<q.error({0,0,0,1})<<endl;

 for (auto v:ppmmesh.vertices()) 
  {
      pmp::Quadric Q=pmp::Quadric();
      Quadric q({0,0,0},{0,0,0});
      for(auto f:ppmmesh.faces(v))
      {
       Q+=pmp::Quadric(normales[f], points[v]);
       q.add(Quadric({normales[f][0],normales[f][1],normales[f][2]},{points[v][0],points[v][1],points[v][1]}));
       
      }
     cout <<Q(points[v])<<" "<<q.error({points[v][0],points[v][1],points[v][1],1})<<endl;
  
  }
 
 
 
     
}




double  compute_diagonalLength(const pmp::SurfaceMesh & ppmmesh, double& xmin ,double& xmax,double& ymin ,double& ymax,double& zmin ,double& zmax)
{
     pmp::VertexProperty<pmp::Point>  points = ppmmesh.get_vertex_property<pmp::Point>("v:point");

     for (auto v : ppmmesh.vertices())     
     {
          auto p= points[v];
          if(p[0]<xmin)
          xmin=p[0];
          if(p[1]<ymin)
          ymin=p[1];
          if(p[2]<zmin)
          zmin=p[2];

          if(p[0]>xmax)
          xmax=p[0];
          if(p[1]>ymax)
          ymax=p[1];
          if(p[2]>zmax)
          zmax=p[2];
          
     }

     return sqrt(pow(xmax-xmin,2)+pow(ymax-ymin,2)+pow(zmax-zmin,2));

}

void computeSaillancy(pmp::SurfaceMesh & ppmmesh)
{
     pmp::Curvature c= pmp::Curvature(ppmmesh);
     c.analyze_tensor();

     
     pmp::VertexProperty<pmp::Point>  points = ppmmesh.get_vertex_property<pmp::Point>("v:point");

     auto saillancyi = ppmmesh.add_vertex_property<vector<pmp::Scalar>>("added:saillancyi");

     auto smoothsaillancy = ppmmesh.add_vertex_property<pmp::Scalar>("added:smoothsaillancy");

     vector<double> max_saillancyi;

     for(int i=0;i<5;i++)
          max_saillancyi.push_back(std::numeric_limits<double>::min());

    
     double xmin ,xmax,ymin,ymax,zmin ,zmax;
     xmin=ymin=zmin=std::numeric_limits<double>::max();
     xmax=ymax=zmax=std::numeric_limits<double>::min();
     double diagonalLength=compute_diagonalLength(ppmmesh,xmin,xmax,ymin,ymax,zmin,zmax);


     double epsilon=2*0.003*diagonalLength;
     vector<int> indicesepsilon={2,3,4,5,6,8,10,12};
     for (auto v : ppmmesh.vertices())     
     {
       
          vector<pmp::Vertex> voisinnages2;
          vector<pmp::Vertex> voisinnages3;
          vector<pmp::Vertex> voisinnages4;
          vector<pmp::Vertex> voisinnages5;
          vector<pmp::Vertex> voisinnages6;
          vector<pmp::Vertex> voisinnages8;
          vector<pmp::Vertex> voisinnages10;
          vector<pmp::Vertex> voisinnages12;
          for(auto vn: ppmmesh.vertices(v))
          {
               
               double distance=sqrt((pow(points[vn][0]-points[v][0],2)+pow(points[vn][1]-points[v][1],2)+pow(points[vn][2]-points[v][2],2)));
               if(distance<2*epsilon)
                    voisinnages2.push_back(vn);
               if(distance<3*epsilon)
                    voisinnages3.push_back(vn);  
               if(distance<4*epsilon)
                    voisinnages4.push_back(vn);
               if(distance<5*epsilon)
                    voisinnages5.push_back(vn); 
               if(distance<6*epsilon)
                    voisinnages6.push_back(vn);
               if(distance<8*epsilon)
                    voisinnages8.push_back(vn);  
               if(distance<10*epsilon)
                    voisinnages10.push_back(vn);
               if(distance<12*epsilon)
                    voisinnages12.push_back(vn); 
          }
          vector<vector<pmp::Vertex>> finalv;
          finalv.push_back(voisinnages2);
          finalv.push_back(voisinnages3);
          finalv.push_back(voisinnages4);
          finalv.push_back(voisinnages5);
          finalv.push_back(voisinnages6);
          finalv.push_back(voisinnages8);
          finalv.push_back(voisinnages10);
          finalv.push_back(voisinnages12);
       
 
          epsilon=0.003*diagonalLength;
          vector<pmp::Scalar> gaussianWeigths;
          for(int i=0;i<8;i++)  
          {
          
          double numerateur=0;   
          double denominateur =0;

          double sigma=pow(indicesepsilon[i]*epsilon,2);
          
          for(int j=0;j<(int)finalv[i].size();j++)
          {
               double norme=0;
               double expo=0;
               norme=(pow(points[finalv[i][j]][0]-points[v][0],2)+pow(points[finalv[i][j]][1]-points[v][1],2)+pow(points[finalv[i][j]][2]-points[v][2],2));
               expo=exp(-norme/2*sigma);
               numerateur+=c.mean_curvature(finalv[i][j])*expo;
               denominateur+=expo;
               
          }
          
          if(denominateur==0||numerateur==0)
               gaussianWeigths.push_back(0.0);
          else
               gaussianWeigths.push_back((numerateur/denominateur));

          }
           vector<pmp::Scalar> saillancy;
           double s1,s2,s3,s4,s5;
                                 
          s1=abs(gaussianWeigths[0]-gaussianWeigths[2]);
          s2=abs(gaussianWeigths[1]-gaussianWeigths[4]);
          s3=abs(gaussianWeigths[2]-gaussianWeigths[5]);
          s4=abs(gaussianWeigths[3]-gaussianWeigths[6]);
          s5=abs(gaussianWeigths[4]-gaussianWeigths[7]);
          
          max_saillancyi[0]=max( max_saillancyi[0],s1);
          max_saillancyi[1]=max( max_saillancyi[1],s2);
          max_saillancyi[2]=max( max_saillancyi[2],s3);
          max_saillancyi[3]=max( max_saillancyi[3],s4);
          max_saillancyi[4]=max( max_saillancyi[4],s5);

           saillancy.push_back(s1);
           saillancy.push_back(s2);
           saillancy.push_back(s3);
           saillancy.push_back(s4);
           saillancy.push_back(s5);
           
       
          saillancyi[v]= saillancy;   
           
          } 


     for (auto v : ppmmesh.vertices())     
     {
     vector<double> local_max_saillancyi;
     for(int i=0;i<5;i++)
          local_max_saillancyi.push_back( saillancyi[v][i]);
     
     for (auto vn : ppmmesh.vertices(v))     
     {
          for(int i=0;i<5;i++){
          local_max_saillancyi[i] = max(local_max_saillancyi[i], (double)saillancyi[vn][i]);
          
          }
     }
     
     double finalsaillancy=0;
     double saliencySum=0;

     for(int i=0;i<5;i++)
     {
          double   factor = pow(max_saillancyi[i]-local_max_saillancyi[i],2);
         
          finalsaillancy += saillancyi[v][i] * factor;
          saliencySum += factor;
     }
     if(finalsaillancy==0||saliencySum==0)
      smoothsaillancy[v] =0;
     else
      smoothsaillancy[v] =finalsaillancy/saliencySum;

     }



     
}








int main(int argc, char const *argv[])
{
     string filename="../object/meshes/armadillo.off";
     //string filename="../object/meshes/sphere.off";
     pmp::SurfaceMesh ppmmesh;
     pmp::read(ppmmesh,filename);
    

   /* for(int i=0;i< faceproperties.size();i++)

          cout<<faceproperties[i]<<endl;

          return true;*/
    
 arma::vec a={0,0,0};
 arma::vec b={0,1,0};
 arma::vec c={0,0,1};
arma::vec ab=b-a;
arma::vec ac=c-a;

arma::vec cro=arma::cross(ab,ac);
arma::vec cros=arma::cross(ac,ab);
//cout <<ab<<endl;
//cout <<ac<<endl;
cout <<cro<<endl;
//cout <<cros<<endl;
      pmp::Quadric Q=pmp::Quadric();
      Quadric q(cro,a);

    //  cout<<q.error({0,0,0,1})<<endl;




   //  computeQuadrix(ppmmesh);
    // computeSaillancy(ppmmesh);
//     auto smoothsaillancy = ppmmesh.get_vertex_property<pmp::Scalar>("added:smoothsaillancy");
   //  for (auto v:ppmmesh.vertices())  
     //     cout<<smoothsaillancy[v]<<endl;

 


    


 
    

   







/*

     for(int i=0;i< faceproperties.size();i++)arma::vec ab=b-a;

          cout<<faceproperties[i]<<endl;

          return true;
     */
     










/*


     auto gauss =ppmmesh.get_vertex_property<vector<pmp::Scalar>>("added:poidgaussian");

          for (auto v : ppmmesh.vertices())  
          {
          vector<pmp::Scalar> poid=gauss[v];
          for(int i=0;i<poid.size();i++)
          {        
                    // cout <<poid[i]<<endl;
          }
          //  cout<<" -------------------------------------"<<endl;
          }

*/








     //pmp::write(ppmmesh,"output.obj");


     return 0;
}

