#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

arma::cube rinvwish(int n, int v, arma::mat S);
arma::colvec vech(arma::mat A);
arma::mat devech(arma::colvec v, int dim);

#endif
