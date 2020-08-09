#include "../include/MeshConnectivity.h"
#include "../include/ElasticShell.h"

#include <mex.h>
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>

#include <iostream>
#include <vector>
#include <algorithm>

template <class SFF>
void helper(
  const MeshConnectivity &mesh, 
  const Eigen::MatrixXd & origPos,
  const Eigen::MatrixXd & curPos,
  const Eigen::VectorXd &thicknesses,
  const int matid)
{
  // initialize default edge DOFs (edge director angles)
  Eigen::VectorXd edgeDOFs;
  SFF::initializeExtraDOFs(edgeDOFs, mesh, curPos);
  // initialize first fundamental forms to those of input mesh
  std::vector<Eigen::Matrix2d> abar;
  ElasticShell<SFF>::firstFundamentalForms(mesh, origPos, abar);
  // initialize second fundamental forms to rest flat
  std::vector<Eigen::Matrix2d> bbar;
  bbar.resize(mesh.nFaces());
  for (int i = 0; i < mesh.nFaces(); i++)
  {
    bbar[i].setZero();
  }
  MaterialModel<SFF> *mat;
  switch (matid)
  {
  case 0:
      mat = new NeoHookeanMaterial<SFF>(lameAlpha, lameBeta);
      break;
  case 1:
      mat = new StVKMaterial<SFF>(lameAlpha, lameBeta);
      break;
  default:
      assert(false);
  }

  {
    int nverts = (int)curPos.rows();
    int nedges = mesh.nEdges();
    int nedgedofs = SFF::numExtraDOFs;
    int freeDOFs = 3 * nverts + nedgedofs * nedges;
    Eigen::VectorXd derivative;
    std::vector<Eigen::Triplet<double> > hessian;
    double energy = ElasticShell<SFF>::elasticEnergy(
      mesh, curPos, edgeDOFs, mat, thicknesses, abars, bbars, &derivative, &hessian);
  }

  delete mat;
}

void mexFunction(
         int          nlhs,
         mxArray      *plhs[],
         int          nrhs,
         const mxArray *prhs[]
         )
{
  //mexPrintf("Compiled at %s on %s\n",__TIME__,__DATE__);
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  Eigen::MatrixXd V;
  Eigen::MatirxXi F;
  double young = 1e9;
  double poisson = 0.3;
  mexErrMsgTxt(nrhs>=2,"nrhs should be >= 2");
  parse_rhs_double(prhs+0,V);
  parse_rhs_index(prhs+1,F);
  // By default use rest pose
  Eigen::MatrixXd U = V;
  int matid = 0;
  int sffid = 0;
  // Paper's thickness in meters
  double thickness = 1e-4;


  // parse all the default arguments


  MeshConnectivity mesh = MeshConnectivity(F);
  const double lameAlpha = young * poisson / (1.0 - poisson * poisson);
  const double lameBeta = young / 2.0 / (1.0 + poisson);
  // uniform thickness for now
  const Eigen::VectorXd thickness  =
    Eigen::VectorXd::Constant(F.rows(),thickness);
  MaterialModel<SFF> *mat;
  switch (matid)
  {
  case 0:
      mat = new NeoHookeanMaterial<SFF>(lameAlpha, lameBeta);
      break;
  case 1:
      mat = new StVKMaterial<SFF>(lameAlpha, lameBeta);
      break;
  default:
      assert(false);
  }
  double E = 0;
  switch (sffid)
  {
  case 0:
      helper<MidedgeAngleTanFormulation>(mesh, V,U, thicknesses, lameAlpha, lameBeta, matid);
      break;
  case 1:
      helper<MidedgeAngleSinFormulation>(mesh, V,U, thicknesses, lameAlpha, lameBeta, matid);
      break;
  case 2:
      helper<MidedgeAverageFormulation>(mesh, V,U, thicknesses, lameAlpha, lameBeta, matid);
      break;
  default:
      assert(false);
  }
  delete mat;



  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}
