#include <RcppEigen.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace std;
using namespace Eigen;


/////////////////////////////////////////
VectorXi extract(VectorXi x, VectorXi ind)
{
  VectorXi out(ind.sum());
  int j=0;
  for(int i=0; i<ind.size(); i++)
  {
    if(ind(i)!=0)
    {
      out(j)=x(i);
      j=j+1;
    }
  }
  return(out);
}

/////////////////////////////////////////
////////////////////////////////////////

MatrixXd get_cv_error_smooth(MatrixXd T_train, MatrixXd T_valid, MatrixXd Y_train, MatrixXd Y_valid, List y_params, List y_penalty_inv)
{
  VectorXd kappa_set=as<VectorXd>(y_params(1));
  MatrixXd B_vals= as<MatrixXd>(y_params(3)), B_vals_weig=as<MatrixXd>(y_params(5));
  int q=kappa_set.size(), ncol=T_train.cols();
  MatrixXd t_train_mtx(T_train.rows(), ncol+1);
  t_train_mtx.col(0)=MatrixXd::Constant(T_train.rows(),1, 1/sqrt(T_train.rows()));
  t_train_mtx.rightCols(ncol)=T_train;
  MatrixXd t_valid_mtx(T_valid.rows(), ncol+1);
  t_valid_mtx.col(0)=MatrixXd::Constant(T_valid.rows(),1, 1/sqrt(T_train.rows()));
  t_valid_mtx.rightCols(ncol)=T_valid;

  MatrixXd coef_w_0=B_vals_weig*Y_train.transpose()*t_train_mtx, coef_w;
  MatrixXd error(q,ncol), V(coef_w_0.cols(), B_vals.cols()), Y_pred;
  for(int k=0; k<q; k++)
  {
    MatrixXd tmp=as<MatrixXd>(as<List>(y_penalty_inv(k))(0));
    VectorXd tmp_1=tmp*coef_w_0.col(0);
    V.row(0)=tmp_1.transpose()*B_vals;
    for(int ncomp=0;ncomp<ncol;ncomp++)
    {
      tmp_1=as<MatrixXd>(as<List>(y_penalty_inv(k))(ncomp+1))*coef_w_0.col(ncomp+1);
      V.row(ncomp+1)=tmp_1.transpose()*B_vals;
      Y_pred=t_valid_mtx.leftCols(ncomp+2)*V.topRows(ncomp+2);
      error(k,ncomp)=(Y_pred-Y_valid).squaredNorm()/Y_valid.cols();
    }
  }
  return error;
}

///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


List calculate_G_and_Y(List t_x,  List X, MatrixXd Y, List x_params, List all_folds)
{
  int n_curves=X.size(), nsample=Y.rows(), K_cv=all_folds.size();
  List B_main_vals=x_params(3), B_inter_vals=x_params(6);
  MatrixXi inter_mat=MatrixXi::Zero(n_curves+n_curves*(n_curves-1)/2,2);
  int k=0;
  for(int i=0; i<n_curves; i++)
  {
    for(int j=i; j<n_curves; j++)
    {
      inter_mat(k,0)=i;
      inter_mat(k,1)=j;
      k=k+1;
    }
  }
  ////////////calculate the list of G and G_inter for the whole data
  List G_mean_main_list(n_curves), G_mean_inter_list(inter_mat.rows()), G_main_list(n_curves),  G_inter_list(inter_mat.rows()), tmp_list(n_curves);

  for(int i=0; i<n_curves; i++)
  {
    MatrixXd tmp_mat=as<MatrixXd>(X(i)).transpose();
    tmp_list(i)=as<MatrixXd>(B_inter_vals(i))*tmp_mat/tmp_mat.rows();
    tmp_mat=as<MatrixXd>(B_main_vals(i))*tmp_mat/tmp_mat.rows();
    G_mean_main_list(i)=tmp_mat.rowwise().mean();
    G_main_list(i)=(tmp_mat.colwise()-tmp_mat.rowwise().mean())/sqrt(nsample);
  }
  for(int i=0; i<inter_mat.rows(); i++)
  {
    MatrixXd tmp_1=as<MatrixXd>(tmp_list(inter_mat(i,0)));
    MatrixXd tmp_2=as<MatrixXd>(tmp_list(inter_mat(i,1)));
    MatrixXd tmp=MatrixXd::Zero(tmp_1.rows()*tmp_2.rows(), tmp_1.cols());
    for(int j=0; j<tmp_1.rows(); j++)
    {
      for(int k=0; k<tmp_2.rows(); k++)
      {
        tmp.row(j*tmp_2.rows()+k)=tmp_1.row(j).cwiseProduct(tmp_2.row(k));
      }
    }
    G_mean_inter_list(i)=tmp.rowwise().mean();
    G_inter_list(i)=(tmp.colwise()-tmp.rowwise().mean())/sqrt(nsample);
  }

  /////////////////////////
  MatrixXd Y_cent=Y.rowwise()-Y.colwise().mean(), Pi=Y_cent*Y_cent.transpose()/nsample/Y.cols();
  Pi=(Pi+Pi.transpose())/2;
  List Pi_train_list(K_cv), Y_train_list(K_cv), Y_valid_list(K_cv);
  for(int i_cv=0; i_cv<K_cv; i_cv++)
  {
    VectorXi omit;
    omit=as<VectorXi>(all_folds(i_cv)).array()-1;
    int nsample_train=nsample-omit.size();
    MatrixXd Y_valid=MatrixXd::Zero(omit.size(), Y.cols());
    MatrixXd Y_train=MatrixXd::Zero(nsample_train, Y.cols());
    VectorXi ind=VectorXi::Zero(nsample);
    for(int i=0; i<omit.size(); i++)
    {
      ind(omit(i))=1;
      Y_valid.row(i)=Y.row(omit(i));
    }
    int i_train=0;
    for(int i=0; i<nsample; i++)
    {
      if(ind(i)==0)
      {
        Y_train.row(i_train)=Y.row(i);
        i_train=i_train+1;
      }
    }
    Y_cent=Y_train.rowwise()-Y_train.colwise().mean();
    MatrixXd tmp_1=Y_cent*Y_cent.transpose()/nsample_train/Y.cols();
    Pi_train_list(i_cv)=(tmp_1+tmp_1.transpose())/2;
    Y_train_list(i_cv)=Y_train;
    Y_valid_list(i_cv)=Y_valid;
  }
  return List::create( _["inter_mat"]=inter_mat, _["G_mean_main_list"]=G_mean_main_list, _["G_mean_inter_list"]=G_mean_inter_list,
                       _["G_main_list"]=G_main_list, _["G_inter_list"]=G_inter_list,  _["Pi"]=Pi, _["Pi_train_list"]=Pi_train_list,
                         _["Y_train_list"]=Y_train_list, _["Y_valid_list"]=Y_valid_list);
}


////////////////////////////////////////////////////////
///////////////////////////////////////////////////////
MatrixXd cal_R_inv(MatrixXd v, List weights, List x_params, int n_main, int n_inter, double tau)
{
  Map<MatrixXd> v_Xd=Map<MatrixXd>(v.data(), v.rows(), v.cols());
  int q=v.cols();
  MatrixXd J0=as<List>(x_params(4))(0), J2=as<List>(x_params(4))(1);
  MatrixXd J00=as<List>(x_params(7))(0), J20=as<List>(x_params(7))(1), J11=as<List>(x_params(7))(2), J02=as<List>(x_params(7))(3);
  MatrixXd out(v.rows(), v.cols());
  Map<MatrixXd> out_Xd=Map<MatrixXd>(out.data(), out.rows(), out.cols());

  int n_main_basis=(as<VectorXi>(x_params(5)))(0), n_inter_basis=(as<VectorXi>(x_params(5)))(1);
  for(int i=0; i<n_main; i++)
  {

    double w0=as<MatrixXd>(weights(0))(i,0), w2=as<MatrixXd>(weights(0))(i,1);
    MatrixXd tmp=w0*J0+tau*w2*J2;
    Map<MatrixXd> tmp_Xd=Map<MatrixXd>(tmp.data(), tmp.rows(), tmp.cols());
    LLT<MatrixXd> lltOf(tmp_Xd);
    out_Xd.block(i*n_main_basis,0, n_main_basis, q)= lltOf.matrixU().solve(v_Xd.block(i*n_main_basis,0, n_main_basis, q));
  }
  for(int i=0; i<n_inter; i++)
  {
    double w00=as<MatrixXd>(weights(1))(i,0), w20=as<MatrixXd>(weights(1))(i,1), w11=as<MatrixXd>(weights(1))(i,2), w02=as<MatrixXd>(weights(1))(i,3);
    MatrixXd tmp=w00*J00+tau*(w20*J20+w11*J11+w02*J02);
    Map<MatrixXd> tmp_Xd=Map<MatrixXd>(tmp.data(), tmp.rows(), tmp.cols());
    LLT<MatrixXd> lltOf(tmp_Xd);
    out_Xd.block(n_main_basis*n_main+i*n_inter_basis, 0, n_inter_basis, q)= lltOf.matrixU().solve(v_Xd.block(n_main_basis*n_main+i*n_inter_basis, 0, n_inter_basis, q));
  }

  return out;
}
/////////////////////////////////////////////////
MatrixXd cal_R_trans_inv(Map<MatrixXd> v, List weights, List x_params, int n_main, int n_inter, double tau)
{
  int q=v.cols();
  MatrixXd J0=as<List>(x_params(4))(0), J2=as<List>(x_params(4))(1);
  MatrixXd J00=as<List>(x_params(7))(0), J20=as<List>(x_params(7))(1), J11=as<List>(x_params(7))(2), J02=as<List>(x_params(7))(3);
  MatrixXd out(v.rows(), v.cols());
  Map<MatrixXd> out_Xd=Map<MatrixXd>(out.data(), out.rows(), out.cols());

  int n_main_basis=(as<VectorXi>(x_params(5)))(0), n_inter_basis=(as<VectorXi>(x_params(5)))(1);
  for(int i=0; i<n_main; i++)
  {

    double w0=as<MatrixXd>(weights(0))(i,0), w2=as<MatrixXd>(weights(0))(i,1);
    MatrixXd tmp=w0*J0+tau*w2*J2;
    Map<MatrixXd> tmp_Xd=Map<MatrixXd>(tmp.data(), tmp.rows(), tmp.cols());
    LLT<MatrixXd> lltOf(tmp_Xd);
    out_Xd.block(i*n_main_basis,0, n_main_basis, q)= lltOf.matrixL().solve(v.block(i*n_main_basis,0, n_main_basis, q));
  }
  for(int i=0; i<n_inter; i++)
  {
    double w00=as<MatrixXd>(weights(1))(i,0), w20=as<MatrixXd>(weights(1))(i,1), w11=as<MatrixXd>(weights(1))(i,2), w02=as<MatrixXd>(weights(1))(i,3);
    MatrixXd tmp=w00*J00+tau*(w20*J20+w11*J11+w02*J02);
    Map<MatrixXd> tmp_Xd=Map<MatrixXd>(tmp.data(), tmp.rows(), tmp.cols());
    LLT<MatrixXd> lltOf(tmp_Xd);
    out_Xd.block(n_main_basis*n_main+i*n_inter_basis, 0, n_inter_basis, q)= lltOf.matrixL().solve(v.block(n_main_basis*n_main+i*n_inter_basis, 0, n_inter_basis, q));
  }

  return out;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////////



void find_orth_basis(Map<MatrixXd> &orthConst_mtx)
{
  int m=orthConst_mtx.cols();
  orthConst_mtx.col(0)=orthConst_mtx.col(0)/((orthConst_mtx.col(0)).norm());
  for(int i=1;i<m;++i)
  {
    MatrixXd tmp=(orthConst_mtx.col(i-1)).transpose()*orthConst_mtx.rightCols(m-i);
    orthConst_mtx.rightCols(m-i)=orthConst_mtx.rightCols(m-i)-orthConst_mtx.col(i-1)*tmp;
    orthConst_mtx.col(i)=orthConst_mtx.col(i)/((orthConst_mtx.col(i)).norm());
  }
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////////


void  power_eig(Map<MatrixXd> Pi, Map<MatrixXd> M, double & max_value, Map<VectorXd> & beta_current, double tol=1e-16)//function used to calculate the first eigenvalue of Pi which is orthogonal to M
{
  int n=Pi.rows(), p=M.rows();
  VectorXd tmp_ini=VectorXd::Constant(p, 1), b_ini=tmp_ini/tmp_ini.norm(), b_old_ini=VectorXd::Zero(p), tmp1_ini=VectorXd::Zero(p);
  Map<VectorXd> tmp=Map<VectorXd>(tmp_ini.data(), tmp_ini.size());
  Map<VectorXd> b=Map<VectorXd>(b_ini.data(), b_ini.size());
  Map<VectorXd> b_old=Map<VectorXd>(b_old_ini.data(), b_old_ini.size());
  Map<VectorXd> tmp1=Map<VectorXd>(tmp1_ini.data(), tmp1_ini.size());
  int count=0;
  double e1=(b-b_old).norm(), e2=(b+b_old).norm();
  double sqrt_tol=sqrt(tol);
  while(count<40 && e1>sqrt_tol && e2>sqrt_tol)
  {
    count++;
    b_old=b;
    tmp=b-M*(M.transpose()*b);
    tmp1.head(n)=Pi*tmp.head(n);
    b=tmp1-M*(M.transpose()*tmp1);
    max_value=b_old.dot(b);
    b=b/b.norm();
    e1=(b-b_old).norm();
    e2=(b+b_old).norm();
  }
  beta_current=b;

}

/////////////////////////////////////////////////
/////////////////////////////////////////////////////

List cal_comp_without_max(Map<MatrixXd> Pi_eig, Map<MatrixXd> M_basis, int upper_comp, double thresh)
{

  int nsample=Pi_eig.cols();
  find_orth_basis(M_basis);
  int p=M_basis.rows();
  VectorXd D_ini=VectorXd::Zero(p), beta_current_ini=VectorXd::Zero(p);
  Map<VectorXd> D=Map<VectorXd>(D_ini.data(), p);
  Map<VectorXd> beta_current=Map<VectorXd>(beta_current_ini.data(), p);
  VectorXd max_value_vec=VectorXd::Zero(upper_comp);
  MatrixXd beta=MatrixXd::Zero(p, upper_comp);
  MatrixXd tmp_matrix;
  int n_comp;

  for(n_comp=0; n_comp<upper_comp; n_comp++)
  {
    double max_value;
    power_eig(Pi_eig, M_basis, max_value, beta_current);
    max_value_vec(n_comp)=max_value;
    beta.col(n_comp)=beta_current;
    if((max_value_vec(n_comp)/max_value_vec.head(n_comp+1).sum()<thresh))
    {break;}
    if(n_comp==upper_comp-1)
    {break;}
    D=VectorXd::Zero(p);
    D.head(nsample)=beta_current.head(nsample);
    D=D/D.norm();
    D=D-M_basis*(M_basis.transpose()*D);
    D=D/D.norm();
    MatrixXd tmp_2=M_basis;
    tmp_matrix=MatrixXd::Zero(M_basis.rows(), M_basis.cols()+1);
    tmp_matrix.leftCols(M_basis.cols())=tmp_2;
    tmp_matrix.col(M_basis.cols())=D;
    new (&M_basis) Map<MatrixXd>(tmp_matrix.data(), tmp_matrix.rows(), tmp_matrix.cols());
  }
  return List::create( _["beta"]=beta.leftCols(n_comp+1),  _["max_value"]=max_value_vec);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////////

List cal_comp_with_max(Map<MatrixXd> Pi_eig, Map<MatrixXd> M_basis, int max_comp)
{
  int nsample=Pi_eig.cols();
  find_orth_basis(M_basis);
  int p=M_basis.rows();
  VectorXd D_ini=VectorXd::Zero(p), beta_current_ini=VectorXd::Zero(p);
  Map<VectorXd> D=Map<VectorXd>(D_ini.data(), p);
  Map<VectorXd> beta_current=Map<VectorXd>(beta_current_ini.data(), p);
  VectorXd max_value_vec=VectorXd::Zero(max_comp);
  MatrixXd beta=MatrixXd::Zero(p, max_comp);
  MatrixXd tmp_matrix;
  int n_comp;
  for(n_comp=0; n_comp<max_comp; n_comp++)
  {
    double max_value;
    power_eig(Pi_eig, M_basis, max_value, beta_current);
    max_value_vec(n_comp)=max_value;
    beta.col(n_comp)=beta_current;
    if(n_comp==max_comp-1)
    {break;}
    D=VectorXd::Zero(p);
    D.head(nsample)=beta_current.head(nsample);
    D=D/D.norm();
    D=D-M_basis*(M_basis.transpose()*D);
    D=D/D.norm();
    MatrixXd tmp_2=M_basis;
    tmp_matrix=MatrixXd::Zero(M_basis.rows(), M_basis.cols()+1);
    tmp_matrix.leftCols(M_basis.cols())=tmp_2;
    tmp_matrix.col(M_basis.cols())=D;
    new (&M_basis) Map<MatrixXd>(tmp_matrix.data(), tmp_matrix.rows(), tmp_matrix.cols());

  }
  return List::create( _["beta"]=beta,  _["max_value"]=max_value_vec);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////////


List cv_over_tau(List list_for_G_Pi, List weight_list, List x_params, List y_penalty_inv, List y_params, List all_folds, VectorXi main_index, VectorXi inter_index, int upper_comp, double thresh)
{

  MatrixXi inter_mat=as<MatrixXi>(list_for_G_Pi(0));
  int n_main=main_index.sum(), n_inter=inter_index.sum();
  int n_curves=x_params(0);
  VectorXi main_all=VectorXi::LinSpaced(n_curves,0,n_curves-1), inter_all=VectorXi::LinSpaced(inter_mat.rows(), 0, inter_mat.rows()-1);
  VectorXd tau_set=as<VectorXd>(x_params(8));
  VectorXi main_effect=extract(main_all, main_index), inter_effect=extract(inter_all, inter_index);
  int n_main_basis=as<VectorXd>(x_params(5))(0), n_inter_basis=as<VectorXd>(x_params(5))(1);
  int nsample=x_params(2);
  int opt_K, opt_tau_index, opt_kappa_index, opt_lambda_index;
  double opt_kappa, opt_lambda, min_error = 1e20, opt_tau;
  MatrixXd opt_z, opt_T;
  VectorXd opt_max_value;
  List opt_errors;
  MatrixXd G_ini(n_main*n_main_basis+n_inter*n_inter_basis, nsample);
  Map<MatrixXd> G=Map<MatrixXd>(G_ini.data(), G_ini.rows(), nsample);


  for(int i=0; i<n_main; i++)
  {
    G.block(i*n_main_basis, 0, n_main_basis, nsample)=as<MatrixXd>(as<List>(list_for_G_Pi["G_main_list"])(main_effect(i)));
  }
  for(int i=0; i<n_inter; i++)
  {
    G.block(n_main_basis*n_main+i*n_inter_basis, 0, n_inter_basis, nsample)=as<MatrixXd>(as<List>(list_for_G_Pi["G_inter_list"])(inter_effect(i)));
  }

  for(int i_tau=0; i_tau<tau_set.size(); i_tau++)
  {
    double tau=tau_set(i_tau);
    MatrixXd R_inv_tran_G=cal_R_trans_inv(G, weight_list, x_params, n_main, n_inter, tau);
    VectorXd lambda_set=as<VectorXd>(x_params(1));
    int  size_lambda=lambda_set.size();
    VectorXi max_comp(size_lambda);
    List max_value_list(size_lambda), z_list(size_lambda), T_list(size_lambda), tmp;

    double lambda, tmp_value;
    MatrixXd Pi=list_for_G_Pi["Pi"];

    for(int i_lambda=0; i_lambda<size_lambda; i_lambda++)
    {
      lambda=lambda_set(i_lambda);
      MatrixXd M(nsample+R_inv_tran_G.rows(), nsample);
      Map<MatrixXd> M_basis=Map<MatrixXd>(M.data(), M.rows(), M.cols());
      M_basis.topRows(nsample)=sqrt(lambda)*(MatrixXd::Identity(nsample, nsample));
      M_basis.bottomRows(R_inv_tran_G.rows())=-R_inv_tran_G;
      Map<MatrixXd> Pi_eig=Map<MatrixXd>(Pi.data(), Pi.rows(), Pi.cols());
      List fit=cal_comp_without_max(Pi_eig, M_basis, upper_comp, thresh);

      MatrixXd beta=as<MatrixXd>(fit["beta"]);
      MatrixXd w=beta.topRows(nsample);
      MatrixXd v=beta.bottomRows(beta.rows()-nsample);
      MatrixXd z=pow(lambda, -0.5)*cal_R_inv(v, weight_list, x_params, n_main, n_inter, tau);
      max_comp(i_lambda)=beta.cols();
      max_value_list(i_lambda)=as<VectorXd>(fit["max_value"]);
      z_list(i_lambda)=z;
      T_list(i_lambda)=w;
    }
    ///conduct cross-validation
    int K_cv=all_folds.size();
    VectorXd kappa_set=as<VectorXd>(y_params(1));
    List  errors(size_lambda);



    for(int fold_ind=0; fold_ind<K_cv; fold_ind++)
    {
      MatrixXd Y_train=as<MatrixXd>(as<List>(list_for_G_Pi["Y_train_list"])(fold_ind));
      MatrixXd Y_valid=as<MatrixXd>(as<List>(list_for_G_Pi["Y_valid_list"])(fold_ind));
      MatrixXd Pi=as<MatrixXd>(as<List>(list_for_G_Pi["Pi_train_list"])(fold_ind));
      VectorXi omit=as<VectorXi>(all_folds(fold_ind)).array()-1;
      int nsample_train=nsample-omit.size();
      MatrixXd R_inv_tran_G_valid(R_inv_tran_G.rows(), omit.size());
      MatrixXd R_inv_tran_G_train(R_inv_tran_G.rows(), nsample_train);
      VectorXi ind=VectorXi::Zero(nsample);
      for(int i=0; i<omit.size(); i++)
      {
        ind(omit(i))=1;
        R_inv_tran_G_valid.col(i)=R_inv_tran_G.col(omit(i));

      }
      int i_train=0;
      for(int i=0; i<nsample; i++)
      {
        if(ind(i)==0)
        {
          R_inv_tran_G_train.col(i_train)=R_inv_tran_G.col(i);
          i_train=i_train+1;
        }
      }
      R_inv_tran_G_train=(R_inv_tran_G_train.colwise()-R_inv_tran_G_train.rowwise().mean())*sqrt(nsample)/sqrt(nsample_train);
      R_inv_tran_G_valid=(R_inv_tran_G_valid.colwise()-R_inv_tran_G_train.rowwise().mean())*sqrt(nsample)/sqrt(nsample_train);

      for(int i_lambda=0; i_lambda<size_lambda; i_lambda++)
      {
        lambda=lambda_set(i_lambda);
        MatrixXd M(nsample_train+R_inv_tran_G_train.rows(), nsample_train);
        Map<MatrixXd> M_basis=Map<MatrixXd>(M.data(), M.rows(), M.cols());
        M_basis.topRows(nsample_train)=sqrt(lambda)*(MatrixXd::Identity(nsample_train, nsample_train));
        M_basis.bottomRows(R_inv_tran_G_train.rows())=-R_inv_tran_G_train;

        Map<MatrixXd> Pi_eig=Map<MatrixXd>(Pi.data(), Pi.rows(), Pi.cols());
        List fit=cal_comp_with_max(Pi_eig, M_basis, max_comp(i_lambda));
        MatrixXd beta=as<MatrixXd>(fit["beta"]);
        MatrixXd T_train=beta.topRows(nsample_train);
        MatrixXd T_valid=pow(lambda, -0.5)*R_inv_tran_G_valid.transpose()*beta.bottomRows(beta.rows()-nsample_train);

        for(int i=0; i<T_train.cols(); i++)
        {
          tmp_value=(T_train.col(i)).norm();
          T_train.col(i)=T_train.col(i)/tmp_value;
          T_valid.col(i)=T_valid.col(i)/tmp_value;
        }

        if(fold_ind==0)
        {
          errors(i_lambda)=get_cv_error_smooth(T_train, T_valid, Y_train, Y_valid, y_params,  y_penalty_inv);
        }else
        {
          errors(i_lambda)=as<MatrixXd>(errors(i_lambda))+get_cv_error_smooth(T_train, T_valid, Y_train, Y_valid, y_params,   y_penalty_inv);
        }
      }
    }


    int tmp_opt_K=0, tmp_opt_lambda_index=0, tmp_opt_kappa_index=0;
    double tmp_opt_kappa=0, tmp_opt_lambda=0, tmp_min_error=1e20;
    for(int i_lambda=0; i_lambda<size_lambda; i_lambda++)
    {
      MatrixXd::Index minRow, minCol;
      double tmp_error=(as<MatrixXd>(errors(i_lambda))).minCoeff(&minRow, &minCol);
      if(tmp_error<tmp_min_error)
      {
        tmp_min_error=tmp_error;
        tmp_opt_lambda_index=i_lambda;
        tmp_opt_kappa_index=minRow;
        tmp_opt_kappa=kappa_set(minRow);
        tmp_opt_lambda=lambda_set(i_lambda);
        tmp_opt_K=minCol+1;
      }
      //Rcout << "tmp_min_error=" <<" " <<tmp_min_error << std::endl;
    }
    if(min_error>tmp_min_error)
    {
      min_error=tmp_min_error;
      opt_lambda_index=tmp_opt_lambda_index;
      opt_tau=tau_set(i_tau);
      opt_tau_index=i_tau;
      opt_kappa_index=tmp_opt_kappa_index;
      opt_K= tmp_opt_K;
      opt_kappa=tmp_opt_kappa;
      opt_lambda=tmp_opt_lambda;
      opt_errors=errors;
      MatrixXd opt_z_1=as<MatrixXd>(z_list(opt_lambda_index));
      opt_z=opt_z_1.leftCols(opt_K);
      MatrixXd opt_T_1=as<MatrixXd>(T_list(opt_lambda_index));
      opt_T=opt_T_1.leftCols(opt_K);
      VectorXd opt_max_value_1=as<VectorXd>(max_value_list(opt_lambda_index));
      opt_max_value=opt_max_value_1.head(opt_K);
    }
    //Rcout << "min_error=" <<" " << min_error << std::endl;
    //Rcout << "########################################" <<std::endl;
  }
  return List::create(_["inter_mat"]=inter_mat, _["main_index"]=main_index, _["inter_index"]=inter_index,   _["min_error"]=min_error, _["opt_K"]=opt_K,   _["opt_lambda"]=opt_lambda,
                      _["opt_lambda_index"]=opt_lambda_index, _["opt_tau"]=opt_tau, _["opt_tau_index"]=opt_tau_index, _["opt_kappa"]=opt_kappa, _["opt_kappa_index"]=opt_kappa_index,
                      _["opt_T"]=opt_T, _["opt_z"]=opt_z, _["opt_max_value"]=opt_max_value, _["G_mean_main_list"]=list_for_G_Pi["G_mean_main_list"], _["G_mean_inter_list"]=list_for_G_Pi["G_mean_inter_list"]);
}


////////////////////////////////////////////////
/////////////////////////////////////////////////////


/////////////////////////////////////////////////////
// [[Rcpp::export]]
List C_cv_fix_effects(List t_x, List X, Eigen::MatrixXd Y, Eigen::VectorXi main_index, Eigen::VectorXi inter_index,   List x_params_raw, List x_params, List y_params, List all_folds, int upper_comp, double thresh)
{
  int n_main=main_index.sum(), n_inter=inter_index.sum();
  List list_for_G_Pi=calculate_G_and_Y(t_x,  X, Y, x_params, all_folds);
  List weight_list(2);
  weight_list(0)=MatrixXd::Constant(n_main, 2, 1);
  weight_list(1)=MatrixXd::Constant(n_inter, 4, 1);
  VectorXd weight_y=MatrixXd::Constant(upper_comp+1,1,1);

  MatrixXd A=as<MatrixXd>(as<List>(y_params(6))(0)), B=as<MatrixXd>(as<List>(y_params(6))(1));
  VectorXd kappa_set=as<VectorXd>(y_params(1));
  List y_penalty_inv(kappa_set.size());
  for(int i=0; i<kappa_set.size(); i++)
  {
    List tmp_list(weight_y.size());
    double kappa=kappa_set(i);
    for(int j=0; j<weight_y.size(); j++)
    {
      tmp_list(j)=(A + kappa*weight_y(j)*B).inverse();
    }
    y_penalty_inv(i)=tmp_list;
  }
  Rcout << "**CV procedure for nonadaptive fitting**"  << std::endl;
  Rcout << "**(used to determine the adaptive constants)**"  << std::endl;

  List fit_cv_fix_effects=cv_over_tau(list_for_G_Pi, weight_list, x_params_raw,  y_penalty_inv, y_params, all_folds, main_index,  inter_index, upper_comp, thresh);
  MatrixXd T=fit_cv_fix_effects["opt_T"];
  MatrixXd B_vals= as<MatrixXd>(y_params(3)), B_vals_weig=as<MatrixXd>(y_params(5)), K=as<MatrixXd>(y_params(4));
  for(int i=0; i<T.cols(); i++)
  {
    T.col(i)=T.col(i)/T.col(i).norm();
  }
  MatrixXd t_train_mtx(T.rows(), T.cols()+1);
  t_train_mtx.col(0)=MatrixXd::Constant(T.rows(),1, 1/sqrt(T.rows()));
  t_train_mtx.rightCols(T.cols())=T;

  MatrixXd coef_w_0=B_vals_weig*Y.transpose()*t_train_mtx, coef_w;
  VectorXd  v(coef_w_0.cols());
  int kappa_index=fit_cv_fix_effects["opt_kappa_index"];
  for(int ncomp=0;ncomp<v.size();ncomp++)
  {
    VectorXd tmp=as<MatrixXd>(as<List>(y_penalty_inv(kappa_index))(ncomp))*coef_w_0.col(ncomp);
    v(ncomp)=tmp.dot(K*tmp);
  }

  weight_y.head(v.size())=v(0)*v.cwiseInverse();
  for(int j=v.size(); j<weight_y.size(); j++)
  {
    weight_y(j)=weight_y(v.size()-1);
  }

  for(int i=0; i<kappa_set.size(); i++)
  {
    List tmp_list(weight_y.size());
    double kappa=kappa_set(i);
    for(int j=0; j<weight_y.size(); j++)
    {
      tmp_list(j)=(A + kappa*weight_y(j)*B).inverse();
    }
    y_penalty_inv(i)=tmp_list;
  }
  MatrixXd J0=as<List>(x_params(4))(0), J2=as<List>(x_params(4))(1);
  MatrixXd J00=as<List>(x_params(7))(0), J20=as<List>(x_params(7))(1), J11=as<List>(x_params(7))(2), J02=as<List>(x_params(7))(3);


  MatrixXd Z=as<MatrixXd>(fit_cv_fix_effects["opt_z"]);
  VectorXd max_value=as<VectorXd>(fit_cv_fix_effects["opt_max_value"]);
  max_value=max_value/max_value.maxCoeff();
  int  n_main_basis=(as<VectorXd>(x_params(5)))(0), n_inter_basis=(as<VectorXd>(x_params(5)))(1);

  if(n_main+n_inter>1)
  {
    MatrixXd tmp=Z.topRows(n_main_basis);
    VectorXd tmp_3=((tmp.transpose()*(J0*tmp)).diagonal()).cwiseAbs();
    double a0=max_value.dot(tmp_3);
    MatrixXd w0=MatrixXd::Constant(n_main, 2,1);
    MatrixXd w1=MatrixXd::Constant(n_inter, 4, 1);
    for(int k=1; k<n_main; k++)
    {
      tmp=Z.block(k*n_main_basis, 0, n_main_basis, Z.cols());
      tmp_3=((tmp.transpose()*(J0*tmp)).diagonal()).cwiseAbs();
      w0(k,0)=a0/(max_value.dot(tmp_3)+a0*1e-6);
    }
    for(int k=0; k<n_inter; k++)
    {
      tmp=Z.block(n_main*n_main_basis+k*n_inter_basis, 0, n_inter_basis, Z.cols());
      tmp_3=((tmp.transpose()*(J00*tmp)).diagonal()).cwiseAbs();
      w1(k,0)=a0/(max_value.dot(tmp_3)+a0*1e-6);
    }
    //////////////////////////////////////
    tmp=Z.topRows(n_main_basis);
    tmp_3=((tmp.transpose()*(J2*tmp)).diagonal()).cwiseAbs();
    a0=max_value.dot(tmp_3);

    for(int k=1; k<n_main; k++)
    {
      tmp=Z.block(k*n_main_basis, 0, n_main_basis, Z.cols());
      tmp_3=((tmp.transpose()*(J2*tmp)).diagonal()).cwiseAbs();
      w0(k,1)=a0/(max_value.dot(tmp_3)+a0*1e-6);
    }
    for(int k=0; k<n_inter; k++)
    {
      tmp=Z.block(n_main*n_main_basis+k*n_inter_basis, 0, n_inter_basis, Z.cols());
      tmp_3=((tmp.transpose()*(J20*tmp)).diagonal()).cwiseAbs();
      w1(k,1)=a0/(max_value.dot(tmp_3)+a0*1e-6);
    }
    for(int k=0; k<n_inter; k++)
    {
      tmp=Z.block(n_main*n_main_basis+k*n_inter_basis, 0, n_inter_basis, Z.cols());
      tmp_3=((tmp.transpose()*(J11*tmp)).diagonal()).cwiseAbs();
      w1(k,2)=a0/(max_value.dot(tmp_3)+a0*1e-6);
    }
    for(int k=0; k<n_inter; k++)
    {
      tmp=Z.block(n_main*n_main_basis+k*n_inter_basis, 0, n_inter_basis, Z.cols());
      tmp_3=((tmp.transpose()*(J02*tmp)).diagonal()).cwiseAbs();
      w1(k,3)=a0/(max_value.dot(tmp_3)+a0*1e-6);
    }

    weight_list(0)=w0;
    weight_list(1)=w1;
  }
  Rcout << "**CV procedure for adaptive fitting**"  << std::endl;

  fit_cv_fix_effects=cv_over_tau(list_for_G_Pi, weight_list, x_params,  y_penalty_inv, y_params, all_folds, main_index,  inter_index, upper_comp, thresh);

  return  List::create(_["fit_cv_fix_effects"]=fit_cv_fix_effects, _["y_penalty_inv"]=y_penalty_inv);

}



/////////////////////////////////////////
/////////////////////////////////////////////////////


double cv_with_fixed_params(List list_for_G_Pi, List fit_cv_fix_effects, List weight_list, List all_folds, List x_params, List y_penalty_inv, List y_params,  VectorXi main_index, VectorXi inter_index)
{
  MatrixXi inter_mat=as<MatrixXi>(list_for_G_Pi(0));
  int n_main=main_index.sum(), n_inter=inter_index.sum();
  int n_curves=x_params(0);
  VectorXi main_all=VectorXi::LinSpaced(n_curves,0,n_curves-1), inter_all=VectorXi::LinSpaced(inter_mat.rows(), 0, inter_mat.rows()-1);
  VectorXd tau_set=as<VectorXd>(x_params(8));
  VectorXi main_effect=extract(main_all, main_index), inter_effect=extract(inter_all, inter_index);
  int n_main_basis=as<VectorXd>(x_params(5))(0), n_inter_basis=as<VectorXd>(x_params(5))(1);


  int nsample=x_params(2);
  ///conduct cross-validation

  int K_cv=all_folds.size();
  VectorXd kappa_set=as<VectorXd>(y_params(1));
  int kappa_index=fit_cv_fix_effects["opt_kappa_index"];

  MatrixXd B_vals= as<MatrixXd>(y_params(3)), B_vals_weig=as<MatrixXd>(y_params(5));
  MatrixXd Y_train, Y_valid, tmp_mat, G_valid, R_inv_tran_G_valid, T_train, T_valid;
  int i_tau=fit_cv_fix_effects["opt_tau_index"];
  int opt_K= fit_cv_fix_effects["opt_K"];
  double error=0;
  MatrixXd G_ini(n_main*n_main_basis+n_inter*n_inter_basis, nsample);
  Map<MatrixXd> G=Map<MatrixXd>(G_ini.data(), G_ini.rows(), nsample);

  for(int i=0; i<n_main; i++)
  {
    G.block(i*n_main_basis, 0, n_main_basis, nsample)=as<MatrixXd>(as<List>(list_for_G_Pi["G_main_list"])(main_effect(i)));
  }
  for(int i=0; i<n_inter; i++)
  {
    G.block(n_main_basis*n_main+i*n_inter_basis, 0, n_inter_basis, nsample)=as<MatrixXd>(as<List>(list_for_G_Pi["G_inter_list"])(inter_effect(i)));
  }

  double tau=tau_set(i_tau);
  MatrixXd R_inv_tran_G=cal_R_trans_inv(G, weight_list, x_params, n_main, n_inter, tau);
  ///calculate the maximum numbers of components
  for(int fold_ind=0; fold_ind<K_cv; fold_ind++)
  {
    MatrixXd Y_train=as<MatrixXd>(as<List>(list_for_G_Pi["Y_train_list"])(fold_ind));
    MatrixXd Y_valid=as<MatrixXd>(as<List>(list_for_G_Pi["Y_valid_list"])(fold_ind));
    MatrixXd Pi=as<MatrixXd>(as<List>(list_for_G_Pi["Pi_train_list"])(fold_ind));

    VectorXi omit=as<VectorXi>(all_folds(fold_ind)).array()-1;
    int nsample_train=nsample-omit.size();
    MatrixXd R_inv_tran_G_valid(R_inv_tran_G.rows(), omit.size());
    MatrixXd R_inv_tran_G_train(R_inv_tran_G.rows(), nsample_train);
    VectorXi ind=VectorXi::Zero(nsample);
    for(int i=0; i<omit.size(); i++)
    {
      ind(omit(i))=1;
      R_inv_tran_G_valid.col(i)=R_inv_tran_G.col(omit(i));

    }
    int i_train=0;
    for(int i=0; i<nsample; i++)
    {
      if(ind(i)==0)
      {
        R_inv_tran_G_train.col(i_train)=R_inv_tran_G.col(i);
        i_train=i_train+1;
      }
    }
    R_inv_tran_G_train=(R_inv_tran_G_train.colwise()-R_inv_tran_G_train.rowwise().mean())*sqrt(nsample)/sqrt(nsample_train);
    R_inv_tran_G_valid=(R_inv_tran_G_valid.colwise()-R_inv_tran_G_train.rowwise().mean())*sqrt(nsample)/sqrt(nsample_train);

    double lambda=fit_cv_fix_effects["opt_lambda"];

    MatrixXd M(nsample_train+R_inv_tran_G_train.rows(), nsample_train);
    Map<MatrixXd> M_basis=Map<MatrixXd>(M.data(), M.rows(), M.cols());
    M_basis.topRows(nsample_train)=sqrt(lambda)*(MatrixXd::Identity(nsample_train, nsample_train));
    M_basis.bottomRows(R_inv_tran_G_train.rows())=-R_inv_tran_G_train;

    Map<MatrixXd> Pi_eig=Map<MatrixXd>(Pi.data(), Pi.rows(), Pi.cols());
    List fit=cal_comp_with_max(Pi_eig, M_basis, opt_K);

    MatrixXd beta=as<MatrixXd>(fit["beta"]);
    MatrixXd T_train=beta.topRows(nsample_train);
    MatrixXd T_valid=pow(lambda, -0.5)*R_inv_tran_G_valid.transpose()*beta.bottomRows(beta.rows()-nsample_train);

    int ncol=T_train.cols();
    for(int i=0; i<ncol; i++)
    {
      double tmp_value=(T_train.col(i)).norm();
      T_train.col(i)=T_train.col(i)/tmp_value;
      T_valid.col(i)=T_valid.col(i)/tmp_value;
    }
    MatrixXd t_train_mtx(T_train.rows(), ncol+1);
    t_train_mtx.col(0)=MatrixXd::Constant(T_train.rows(),1, 1/sqrt(T_train.rows()));
    t_train_mtx.rightCols(ncol)=T_train;
    MatrixXd t_valid_mtx(T_valid.rows(), ncol+1);
    t_valid_mtx.col(0)=MatrixXd::Constant(T_valid.rows(),1, 1/sqrt(T_train.rows()));
    t_valid_mtx.rightCols(ncol)=T_valid;


    MatrixXd coef_w_0=B_vals_weig*Y_train.transpose()*t_train_mtx, coef_w;
    MatrixXd  V(coef_w_0.cols(), B_vals.cols()), Y_pred;


    for(int ncomp=0;ncomp<V.rows();ncomp++)
    {
      VectorXd tmp=as<MatrixXd>(as<List>(y_penalty_inv(kappa_index))(ncomp))*coef_w_0.col(ncomp);
      V.row(ncomp)=tmp.transpose()*B_vals;
    }

    Y_pred=t_valid_mtx*V;
    error=error+(Y_pred-Y_valid).squaredNorm()/Y_valid.cols();
  }
  return error;
}


////////////////////////////////////////
// [[Rcpp::export]]
Eigen::MatrixXd C_pred_ff_inter(List fit_cv, Eigen::MatrixXd Y_train, List X_test, List x_params, List y_params, List y_penalty_inv)
{
  MatrixXi inter_mat=fit_cv["inter_mat"];
  int kappa_index=fit_cv["opt_kappa_index"];
  MatrixXd z=fit_cv["opt_z"], T=fit_cv["opt_T"];
  VectorXi main_index=fit_cv["main_index"], inter_index=fit_cv["inter_index"];
  int n_main=main_index.sum(), n_inter=inter_index.sum();
  int n_curves=x_params(0), n_test=(as<MatrixXd>(X_test(0))).rows();
  VectorXi main_all=VectorXi::LinSpaced(n_curves,0,n_curves-1), inter_all=VectorXi::LinSpaced(inter_mat.rows(), 0, inter_mat.rows()-1);
  VectorXd tau_set=as<VectorXd>(x_params(8));
  VectorXi main_effect=extract(main_all, main_index), inter_effect=extract(inter_all, inter_index);


  int n_main_basis=as<VectorXi>(x_params(5))(0), n_inter_basis=as<VectorXi>(x_params(5))(1);
  int nsample=x_params(2);
  List B_main_vals=x_params(3), B_inter_vals=x_params(6);
  MatrixXd G_test(n_main*n_main_basis+n_inter*n_inter_basis, n_test);
  MatrixXd B_vals= as<MatrixXd>(y_params(3)), B_vals_weig=as<MatrixXd>(y_params(5));
  for(int i=0; i<n_main; i++)
  {
    MatrixXd tmp_mat=as<MatrixXd>(X_test(main_effect(i))).transpose();
    tmp_mat=as<MatrixXd>(B_main_vals(main_effect(i)))*tmp_mat/tmp_mat.rows();
    VectorXd tmp_2=as<VectorXd>(as<List>(fit_cv["G_mean_main_list"])(main_effect(i)));
    G_test.block(i*n_main_basis, 0, n_main_basis, n_test)= (tmp_mat.colwise()-tmp_2)/sqrt(nsample);
  }
  for(int i=0; i<n_inter; i++)
  {
    MatrixXd tmp_1=as<MatrixXd>(X_test(inter_mat(inter_effect(i),0))).transpose();
    tmp_1=as<MatrixXd>(B_inter_vals(inter_mat(inter_effect(i),0)))*tmp_1/tmp_1.rows();
    MatrixXd tmp_2=as<MatrixXd>(X_test(inter_mat(inter_effect(i),1))).transpose();
    tmp_2=as<MatrixXd>(B_inter_vals(inter_mat(inter_effect(i),1)))*tmp_2/tmp_2.rows();
    MatrixXd tmp(tmp_1.rows()*tmp_2.rows(), tmp_1.cols());
    for(int j=0; j<tmp_1.rows(); j++)
    {
      for(int k=0; k<tmp_2.rows(); k++)
      {
        tmp.row(j*tmp_2.rows()+k)=tmp_1.row(j).cwiseProduct(tmp_2.row(k));
      }
    }
    VectorXd tmp_3=as<VectorXd>(as<List>(fit_cv["G_mean_inter_list"])(inter_effect(i)));
    G_test.block(n_main_basis*n_main+i*n_inter_basis, 0, n_inter_basis, n_test)= (tmp.colwise()-tmp_3)/sqrt(nsample);
  }
  MatrixXd T_test=G_test.transpose()*z;
  for(int i=0; i<T.cols(); i++)
  {
    double tmp_value=T.col(i).norm();
    T.col(i)=T.col(i)/tmp_value;
    T_test.col(i)=T_test.col(i)/tmp_value;
  }
  int ncol=T.cols();
  MatrixXd t_train_mtx(T.rows(), ncol+1);
  t_train_mtx.col(0)=MatrixXd::Constant(T.rows(),1, 1/sqrt(T.rows()));
  t_train_mtx.rightCols(ncol)=T;
  MatrixXd t_test_mtx(T_test.rows(), ncol+1);
  t_test_mtx.col(0)=MatrixXd::Constant(T_test.rows(),1, 1/sqrt(T.rows()));
  t_test_mtx.rightCols(ncol)=T_test;
  MatrixXd coef_w_0=B_vals_weig*Y_train.transpose()*t_train_mtx, coef_w;
  MatrixXd  V(coef_w_0.cols(), B_vals.cols()), Y_pred;


  for(int ncomp=0;ncomp<V.rows();ncomp++)
  {
    MatrixXd tmp=as<MatrixXd>(as<List>(y_penalty_inv(kappa_index))(ncomp))*coef_w_0.col(ncomp);
    V.row(ncomp)=tmp.transpose()*B_vals;
  }

  Y_pred=t_test_mtx*V;
  return Y_pred;
}
///////////////////////////////////////////////////////////////////////
//////////////////////



/////////////////////////////////////////////////////
// [[Rcpp::export]]
List C_stepwise_adaptive(List t_x,  List X, Eigen::MatrixXd Y, List x_params_raw, List x_params, List y_params, List all_folds, int upper_comp, double thresh)
{

  int n_curves=X.size(),  n_all_main=n_curves;
  List list_for_G_Pi=calculate_G_and_Y(t_x,  X, Y, x_params_raw, all_folds);
  MatrixXi inter_mat=as<MatrixXi>(list_for_G_Pi(0));
  int n_all_inter=inter_mat.rows();
  VectorXd tau_set=as<VectorXd>(x_params(8));
  VectorXi main_all=VectorXi::LinSpaced(n_curves,0,n_curves-1), inter_all=VectorXi::LinSpaced(inter_mat.rows(), 0, inter_mat.rows()-1);

  ///////start the cv procedure


  int total_n_models=20;

  MatrixXi main_index_saved=MatrixXi::Zero(n_all_main, total_n_models), inter_index_saved=MatrixXi::Zero(n_all_inter,total_n_models);
  VectorXd cv_error_saved(total_n_models);

  double min_cv_error=1e20;
  List fit_cv_fix_effects, fit_opt_effects;
  int count=0;
  VectorXi main_index=VectorXi::Constant(n_all_main,1), inter_index=VectorXi::Zero(n_all_inter);
  VectorXi opt_main_index, opt_inter_index;
  VectorXd  SSR_1=VectorXd::Zero(n_all_main), SSR_2=VectorXd::Zero(n_all_inter);

  VectorXd weight_y=VectorXd::Constant(upper_comp+1, 1);

  MatrixXd A=as<MatrixXd>(as<List>(y_params(6))(0)), B=as<MatrixXd>(as<List>(y_params(6))(1));
  VectorXd kappa_set=as<VectorXd>(y_params(1));
  List y_penalty_inv(kappa_set.size());
  for(int i=0; i<kappa_set.size(); i++)
  {
    List tmp_list(weight_y.size());
    double kappa=kappa_set(i);
    for(int j=0; j<weight_y.size(); j++)
    {
      tmp_list(j)=(A + kappa*weight_y(j)*B).inverse();
    }
    y_penalty_inv(i)=tmp_list;
  }
  while(count<total_n_models)
  {
    Rcout << "################################################################## "  << std::endl;
    Rcout << "Step " <<" " <<  count+1 << std::endl;
    Rcout << "** CV procedure for calculation of CV error** "  << std::endl;

    int n_main=main_index.sum(), n_inter=inter_index.sum();

    List weight_list(2);
    weight_list(0)=MatrixXd::Constant(n_main, 2, 1);
    weight_list(1)=MatrixXd::Constant(n_inter, 4, 1);

    fit_cv_fix_effects=cv_over_tau(list_for_G_Pi, weight_list, x_params_raw,  y_penalty_inv, y_params, all_folds, main_index,  inter_index, upper_comp, thresh);

    double current_cv_error=fit_cv_fix_effects["min_error"];
    main_index_saved.col(count)=main_index;
    inter_index_saved.col(count)=inter_index;
    cv_error_saved(count)=current_cv_error;
    Rcout << "min_cv_error" <<" " <<  min_cv_error << std::endl;
    Rcout << "current_cv_error" <<" " <<  current_cv_error << std::endl;
    if(min_cv_error>current_cv_error)
    {
      fit_opt_effects=fit_cv_fix_effects;
      min_cv_error=current_cv_error;
      opt_main_index=main_index;
      opt_inter_index=inter_index;
    }
    else
    {
      break;
    }
    Rcout << "** determine effects added or removed from the current model ** "  << std::endl;
    //determine next model

    VectorXi tmp_main_index=main_index;
    for(int i_main=0; i_main<n_all_main; i_main++)
    {
      tmp_main_index(i_main)=1-main_index(i_main);
      int ind=0;
      if((tmp_main_index.sum()==0)&&(inter_index.sum()==0))
      {
        ind=1;
        SSR_1(i_main)=1e20;
      }
      for(int j=0;j<count+1;j++)
      {
        if(((tmp_main_index-main_index_saved).col(j).cwiseAbs().sum()==0)&& ((inter_index-inter_index_saved).col(j).cwiseAbs().sum()==0))
        {
          ind=1;
          SSR_1(i_main)=1e20;
          break;
        }
      }
      for(int j=0;j<n_all_inter;j++)
      {
        if(inter_index(j)==1)
        {
          if((tmp_main_index(inter_mat(j,0))==0)||(tmp_main_index(inter_mat(j,1))==0))
          {
            ind=1;
            SSR_1(i_main)=1e20;
            break;
          }
        }
      }
      if(ind==0)
      {
        weight_list(0)=MatrixXd::Constant(tmp_main_index.sum(), 2, 1);
        weight_list(1)=MatrixXd::Constant(inter_index.sum(), 4, 1);
        SSR_1(i_main)=cv_with_fixed_params(list_for_G_Pi, fit_cv_fix_effects, weight_list,  all_folds, x_params_raw, y_penalty_inv, y_params,  tmp_main_index,  inter_index);
      }
      tmp_main_index=main_index;
    }

    VectorXi tmp_inter_index=inter_index;
    for(int i_inter=0; i_inter<n_all_inter; i_inter++)
    {
      tmp_inter_index(i_inter)=1-inter_index(i_inter);
      int ind=0;
      if((main_index.sum()==0)&&(tmp_inter_index.sum()==0))
      {
        ind=1;
        SSR_2(i_inter)=1e20;
      }
      for(int j=0;j<count+1;j++)
      {
        if(((main_index-main_index_saved).col(j).cwiseAbs().sum()==0)&& ((tmp_inter_index-inter_index_saved).col(j).cwiseAbs().sum()==0))
        {
          ind=1;
          SSR_2(i_inter)=1e20;
          //break;
        }
      }
      for(int j=0;j<n_all_inter;j++)
      {
        if(tmp_inter_index(j)==1)
        {
          if((main_index(inter_mat(j,0))==0)||(main_index(inter_mat(j,1))==0))
          {
            ind=1;
            SSR_2(i_inter)=1e20;
            // break;
          }
        }
      }
      if(ind==0)
      {
        weight_list(0)=MatrixXd::Constant(main_index.sum(), 2, 1);
        weight_list(1)=MatrixXd::Constant(tmp_inter_index.sum(), 4, 1);
        SSR_2(i_inter)=cv_with_fixed_params(list_for_G_Pi, fit_cv_fix_effects, weight_list,  all_folds, x_params_raw, y_penalty_inv, y_params,  main_index, tmp_inter_index);
      }
      tmp_inter_index=inter_index;
    }

    VectorXd::Index min_Ind_1, min_Ind_2;
    double min_value_1=SSR_1.minCoeff(&min_Ind_1), min_value_2=SSR_2.minCoeff(&min_Ind_2);
    if((min_value_1>9e19)&&(min_value_2>9e19))
    {
      break;
    }

    if(min_value_1<min_value_2)
    {
      main_index(min_Ind_1)=1-main_index(min_Ind_1);
    }
    else
    {
      inter_index(min_Ind_2)=1-inter_index(min_Ind_2);
    }
    count=count+1;
    if(main_index.sum()>0)
    {
      VectorXi selected_main(main_index.sum());
      int i_select=0;
      for(int i=0; i<main_all.size();i++)
      {
        if(main_index(i)==1)
        {
          selected_main(i_select)=main_all(i)+1;
          i_select=i_select+1;
        }
      }
      Rcout << "the main effecs after this step=" <<" " <<  selected_main.transpose() << std::endl;
    }
    else
    {
      Rcout << "the main effecs after this step is empty!"  << std::endl;
    }
    if(inter_index.sum()>0)
    {
      VectorXi selected_inter(inter_index.sum());
      int i_select=0;
      for(int i=0; i<inter_all.size();i++)
      {
        if(inter_index(i)==1)
        {
          selected_inter(i_select)=inter_all(i);
          i_select=i_select+1;
        }
      }
      Rcout << "the interaction effecs after this step==" << std::endl;
      for(int j=0;j<selected_inter.size(); j++)
      {
        Rcout<< "("<<inter_mat(selected_inter(j),0)+1<<","<<inter_mat(selected_inter(j),1)+1<<")"<<std::endl;
      }
    }
    else
    {
      Rcout << "the interaction effecs after this step is empty!" << std::endl;
    }
  }
  Rcout << "################################################################## "  << std::endl;
  if(opt_main_index.sum()>0)
  {
    VectorXi selected_main(opt_main_index.sum());
    int i_select=0;
    for(int i=0; i<main_all.size();i++)
    {
      if(opt_main_index(i)==1)
      {
        selected_main(i_select)=main_all(i)+1;
        i_select=i_select+1;
      }
    }
    Rcout << "finally selected main effecs=" <<" " << selected_main.transpose()  << std::endl;
  }
  else
  {
    Rcout << "finally selected main effecs is empty"   << std::endl;
  }
  if(opt_inter_index.sum()>0)
  {
    VectorXi selected_inter(opt_inter_index.sum());
    int i_select=0;
    for(int i=0; i<inter_all.size();i++)
    {
      if(opt_inter_index(i)==1)
      {
        selected_inter(i_select)=inter_all(i);
        i_select=i_select+1;
      }
    }
    Rcout << "finally selected interation effecs=" << std::endl;
    for(int j=0;j<selected_inter.size(); j++)
    {
      Rcout<< "("<<inter_mat(selected_inter(j),0)+1<<","<<inter_mat(selected_inter(j),1)+1<<")"<<std::endl;
    }
  }
  else
  {
    Rcout << "finally selected interation effecs is empty!" << std::endl;
  }
  Rcout << "################################################################## "  << std::endl;




  return  List::create(_["opt_main_index"]=opt_main_index, _["opt_inter_index"]=opt_inter_index, _["inter_mat"]=inter_mat);

}

/////////////////////////////////////////
/////////////////////////////////////////
////////////////////////////////////////
// [[Rcpp::export]]
List C_find_coef_ff_interaction(List fit_cv, List X_train, Eigen::MatrixXd Y_train, List x_params, List y_params, List y_penalty_inv)
{
  MatrixXi inter_mat=fit_cv["inter_mat"];
  int kappa_index=fit_cv["opt_kappa_index"];
  MatrixXd z=fit_cv["opt_z"], T=fit_cv["opt_T"];
  VectorXi main_index=fit_cv["main_index"], inter_index=fit_cv["inter_index"];
  int n_main=main_index.sum(), n_inter=inter_index.sum();
  int n_curves=x_params(0);
  VectorXi main_all=VectorXi::LinSpaced(n_curves,0,n_curves-1), inter_all=VectorXi::LinSpaced(inter_mat.rows(), 0, inter_mat.rows()-1);
  VectorXd tau_set=as<VectorXd>(x_params(8));
  VectorXi main_effect=extract(main_all, main_index), inter_effect=extract(inter_all, inter_index);
  int n_main_basis=as<VectorXi>(x_params(5))(0), n_inter_basis=as<VectorXi>(x_params(5))(1);
  int nsample=x_params(2);
  List B_main_vals=x_params(3), B_inter_vals=x_params(6);
  MatrixXd B_vals= as<MatrixXd>(y_params(3)), B_vals_weig=as<MatrixXd>(y_params(5));
  int ncol=T.cols();
  MatrixXd scale_factor=MatrixXd::Identity(ncol, ncol);
  for(int i=0; i<ncol; i++)
  {
    scale_factor(i,i)=1/T.col(i).norm();
    T.col(i)=T.col(i)*scale_factor(i,i);
  }
  MatrixXd t_train_mtx(T.rows(), ncol+1);
  t_train_mtx.col(0)=MatrixXd::Constant(T.rows(),1, 1/sqrt(T.rows()));
  t_train_mtx.rightCols(ncol)=T;
  MatrixXd coef_w_0=B_vals_weig*Y_train.transpose()*t_train_mtx, coef_w;
  MatrixXd  V(coef_w_0.cols(), B_vals.cols());
  for(int ncomp=0;ncomp<V.rows();ncomp++)
  {
    MatrixXd tmp=as<MatrixXd>(as<List>(y_penalty_inv(kappa_index))(ncomp))*coef_w_0.col(ncomp);
    V.row(ncomp)=tmp.transpose()*B_vals;
  }
  VectorXd intercept=(V.row(0)/sqrt(T.rows()));
  MatrixXd W=z*scale_factor*V.bottomRows(V.rows()-1)/sqrt(nsample);
  List coef_main(n_main);
  for(int i=0; i<n_main; i++)
  {
    MatrixXd tmp_1=(as<MatrixXd>(B_main_vals(main_effect(i)))).transpose()*W.block(i*n_main_basis, 0, n_main_basis, W.cols());
    coef_main(i)=tmp_1;
    intercept=intercept-tmp_1.transpose()*(as<MatrixXd>(X_train(main_effect(i)))).colwise().mean().transpose()/tmp_1.rows();
  }
  List coef_inter(n_inter);
  for(int i=0; i<n_inter; i++)
  {
    MatrixXd tmp_1=as<MatrixXd>(B_inter_vals(inter_mat(inter_effect(i),0)));
    MatrixXd tmp_2=as<MatrixXd>(B_inter_vals(inter_mat(inter_effect(i),1)));
    List inter_effects_list(W.cols());
    MatrixXd W_tmp=W.block(n_main_basis*n_main+i*n_inter_basis, 0, n_inter_basis, W.cols());
    for(int t=0; t<W.cols(); t++)
    {
      VectorXd W_row_vec=W_tmp.col(t);
      Map<MatrixXd> tmp_W(W_row_vec.data(), tmp_2.rows(), tmp_1.rows());
      MatrixXd tmp_coef=tmp_1.transpose()*tmp_W.transpose()*tmp_2;
      inter_effects_list(t)=tmp_coef;
      MatrixXd tmp_3=(as<MatrixXd>(X_train(inter_mat(inter_effect(i),0)))*tmp_coef).cwiseProduct(as<MatrixXd>(X_train(inter_mat(inter_effect(i),1))))/tmp_coef.rows()/tmp_coef.cols();
      intercept(t)=intercept(t)-tmp_3.rowwise().sum().mean();
    }

    coef_inter(i)=inter_effects_list;
  }

  MatrixXi inter_effects(inter_effect.size(),inter_mat.cols());
  for(int i=0; i<inter_effects.rows(); i++)
  {
    inter_effects.row(i)=inter_mat.row(inter_effect(i));
  }

  return  List::create(_["intercept"]=intercept, _["coef_main"]=coef_main, _["coef_inter"]=coef_inter,
                       _["main_effects"]=main_effect, _["inter_effects"]=inter_effects);

}

///////////////////////////////////////////////////////////////////////
//////////////////////
