\name{NEWS}
\title{News for Package 'FDboost'}

\section{Changes in FDboost version 0.3-1 (2018-05-10)}{
  \subsection{Bug-fixes}{
    \itemize{
    \item \code{cvrisk} has per default adequate folds for a noncyclic fitted FDboostLSS model,  
          see isssue \href{https://github.com/boost-R/FDboost/issues/14}{#14} 
    }
    }
  \subsection{Miscellaneous}{
    \itemize{
    \item replace cBind which is deprecated with cbind
    } 
  }
}

\section{Changes in FDboost version 0.3-0 (2017-05-31)}{
    \subsection{User-visible changes}{
    \itemize{
    \item new function \code{bootstrapCI()} to compute bootstrapped coefficients 
    \item add the dataset 'emotion' containing EEG and EMG measures under different experimental conditions
    \item with scalar response, \code{FDboost()} works with the response as 
          vector and not as matrix with one row; 
          thus, \code{fitted()} and \code{predict()} return a vector 
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
    \item \code{update.FDboost()} works now with scalar response 
    \item \code{FDboost()} works with family \code{Binomial(type = "glm")}, 
          see isssue \href{https://github.com/boost-R/FDboost/issues/1}{#1} 
    \item \code{applyFolds()} works for factor response, 
          see issue \href{https://github.com/boost-R/FDboost/issues/7}{#7}
    \item \code{cvLong} and \code{cvMA} return a matrix for only one resampling 
           fold with \code{B = 1} (proposed by Almond Stoecker)
    }
    }
   \subsection{Miscellaneous}{
    \itemize{
	  \item adapt \pkg{FDboost} to \pkg{mboost} 2.8-0 that allows for mstop = 0
	  \item restructure FDboostLSS() such that it calls mboostLSS_fit() from \pkg{gamboostLSS} 2.0-0  
	  \item in \pkg{FDboost}, set \code{options("mboost_indexmin" = +Inf)} to disable the 
	        internal use of ties in model fitting, as this breaks some methods for models with response 
	        in long format and for models containing \code{bhistx}, 
	        see issue \href{https://github.com/boost-R/FDboost/issues/10}{#10}
	  \item deprecate \code{validateFDboost()}, 
	        use \code{applyFolds()} and \code{bootstrapCI()} instead
    }
  }
}


\section{Changes in FDboost version 0.2-0 (2016-05-26)}{
    \subsection{User-visible changes}{
    \itemize{
    \item add function applyFolds() to compute the optimal stopping iteration
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
    \item allow for extrapolation in predict() with bbsc()
    }
    }
}

\section{Changes in FDboost version 0.1-2 (2016-04-22)}{
    \subsection{Bug-fixes}{
    \itemize{
    \item bugfix in bolsc(): correctly use index in bolsc() / bbsc(), 
          before: for the computation of Z each observation was used only once 
    }
    }
    \subsection{User-visible changes}{
    \itemize{
    \item add function \%Xa0\% that computes a row-tensor product of two base-learners where 
          the penalty in one direction is zero  
    \item add function reweightData() that computes the data for Bootstrap or cross-falidation folds
	  \item add function stabsel.FDboost() that refits the smooth offset in each fold  
	  \item add argument 'fun' to validateFDboost() 
	  \item add update.FDboost() that overwrites update.mboost()
    }
  }
   \subsection{Miscellaneous}{
    \itemize{
	  \item FDboost() works with family = Binomial()  		 
    }
  }
}


\section{Changes in FDboost version 0.1-1 (2016-04-06)}{
  \subsection{Bug-fixes}{
    \itemize{
	    \item fix oobpred in validateFDboost() for irregular response and resampling on the level of curves 
	    and thus plot.validateFDboost() works for that case 
	    \item fix scope of formula in FDboost(): now the formula given to mboost() within FDboost() uses the variables in the environment of the formula specified in FDboost() 
    }
  }
    \subsection{Miscellaneous}{
    \itemize{
	  \item plot.FDboost() works for more effects, especially for effects like bolsc() \%X\% bhistx()  		 
    }
  }
}


\section{Changes in FDboost version 0.1-0 (2016-03-10)}{
  \subsection{User-visible changes}{
    \itemize{
	    \item new operator \%A0\% for Kronecker product of two base-learners with 
		      anisotropic penalty for the special case where lambda1 or lambda2 is zero 
	    \item the base-learner bbsc() can be used with center = TRUE, derived by Almond Stoecker 
		\item in FDboostLSS() a list of one-sided formulas can be specified for timeformula  	
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item FDboostLSS works with families = GammaLSS()
    }
  }
  \subsection{Miscellaneous}{
    \itemize{
	  \item operator \%A\% uses weights in model call; only works correctly for weights on level 
            of blg1 and blg2 (which is the same as weights on rows and columns of the response matrix)	  
      \item call to internal functions of mboost is done using mboost_intern() 	
      \item hyper_olsc() is based on hyper_ols() of mboost 		 
    }
  }
}

\section{Changes in FDboost version 0.0-17 (2016-02-25)}{
  \subsection{User-visible changes}{
    \itemize{
	    \item changed the operator \%Xc\% for row tensor product of two scalar covariates. 
		  The design matrix of the interaction effects is constrained such that the interaction is 
		  centred around the intercept and around the two main effects of the scalar covariates (experimental!);  
		  use e.g. as bols(x1) \%Xc\% bols(x2)
    }
  }
}

\section{Changes in FDboost version 0.0-16 (2016-02-22)}{
  \subsection{User-visible changes}{
    \itemize{
	    \item changed the operator \%Xc\% for row tensor product where the sum-to-zero constraint is applied to
		  the design matrix resulting from the row-tensor product (experimental!), 
		  such that first a, intercept-column is added to the design-matrix and then the sum-to-zero constraint 
		  is applied, use e.g. as bolsc(x1) \%Xc\% bolsc(x2)
		  \item use the functional index s as argsvals in the FPCA conducted within bfpc()
    }
  }
}

\section{Changes in FDboost version 0.0-15 (2016-02-12)}{
  \subsection{User-visible changes}{
    \itemize{
	    \item new operator \%A\% that implies anisotropic penalties for differently specified df in the two base-learners
    }
  }
	\subsection{Bug-fixes}{
    \itemize{
      \item do not penalize in direction of ONEx in smooth intercept specified implicitly by ~1, as bols(ONEx, intercept=FALSE, df=1) \%A\% bbs(time)	  
    }
  }
   \subsection{Miscellaneous}{
    \itemize{
      \item do not expand an effect that contains \%A\% or \%O\% with the timeformula, allowing for different effects over time for the 
        effects in the model 	  
    }
  }
}

\section{Changes in FDboost version 0.0-14 (2016-02-11)}{
  \subsection{User-visible changes}{
    \itemize{
	    \item add the function FDboostLSS() to fit GAMLSS models with functional data 
		      using R-package gamboostLSS
	    \item new operator \%Xc\% for row tensor product where the sum-to-zero constraint is applied to
		  the design matrix resulting from the row-tensor product (experimental!)
	    \item allow newdata to be a list in predict.FDboost() in combination with signal base-learners 
        \item expand coef.FDboost() such that it works for 3-dimensional tensor products
          of with bhistx() the form bhistx() \%X\% bolsc() \%X\% bolsc() (with David Ruegamer)
        \item add a new possibility for scalar-on-function regression: 
          for timeformula=NULL, no Kronecker-product with 1 is used, which 
          changes the penalty as otherwise in the direction of 1 is penalized as well. 
    }
  }
  \subsection{Miscellaneous}{
    \itemize{
      \item new dependency on R-package gamboostLSS 
	    \item remove dependency on R-package MASS 
      \item use the argument 'prediction' in the internal computation 
         of the base-learners (work in progress) 
      \item throw an error if 'timeLab' of the hmatrix-object in bhistx() is not
         equal to the time-variable in 'timeformula'. 
    }
  }
}


\section{Changes in FDboost version 0.0-13 (2015-11-17)}{
  \subsection{User-visible changes}{
    \itemize{
	\item in function FDboost() the offset is supplied differently, for a scalar offset, use offset = "scalar", the default is still the same offset=NULL
	\item predict.FDboost() has new argument toFDboost (logical)
	\item fitted.FDboost() has argument toFDboost explicitly and not only in ...
    \item new base-learner bhistx() especially suited for effects with \%X\%, like bhistx \%X\% bolsc 
    \item coef.FDboost() and plot.FDboost() suited for effects like bhistx \%X\% bolsc
    \item for predict.FDboost() with effects bhistx() and newdata the latest mboostPatch is necessary
    }
  }
    \subsection{Bug-fixes}{
    \itemize{
      \item check for necessity of smooth offset works for missings in regular response (spotted by Tore Erdmann) 
    }
  }
}

\section{Changes in FDboost version 0.0-12 (2015-09-15)}{
    \itemize{
    \item Internal experimental version.
    } 
}

\section{Changes in FDboost version 0.0-11 (2015-06-01)}{
  \subsection{User-visible changes}{
    \itemize{
    \item integrationWeights() gives equal weights for regular grids
    \item new base-learner bfpc() for a functional covariate where
    functional covariate and the coeffcient are both expanded using fPCA (experimental feature!); 
    only works for regularly observed functional covariate. 
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item the function coef.FDboost() only works for bhist() if the time variable is the same in the timeformula and in bhist()
      \item predict.FDboost() has a check that for newdata only type="link" can be predicted
    }
  }
}

\section{Changes in FDboost version 0.0-10 (2015-04-16)}{
  \subsection{User-visible changes}{
    \itemize{
    \item change the default in difference-penalties to first order difference penalty
    differences=1, as then the effects are better identifiable 
      \item new method cvrisk.FDboost() that uses per default 
      sampling on the levels of curves, which is important for functional response
      \item reorganize documentation of cvrisk() and validateFDboost()
	  \item in bhist(): effect can be standardized 
    }
  }
  \subsection{Miscellaneous}{
    \itemize{
      \item add a CITATION file 
	  \item use mboost 2.4-2 as it exports all important functions 
    }
  }
  \subsection{Bug-fixes}{
    \itemize{
      \item main argument is always passed in plot.FDboost()
      \item bhist() and bconcurrent() work for equal time and s
      \item predict.FDboost() works with tensor-product base-learners bl1 \%X\% bl2
    }
  }
}


