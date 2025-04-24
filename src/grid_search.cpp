#include <Rcpp.h>
using namespace Rcpp;

// All the h/g-functions

double finv(double f) {
  return f / sqrt(1 + pow(f,2));
}

double f(double r) {
  return r / sqrt(1 - pow(r,2));
}

double g2(double p1, double p3, double c1) {
  return (p3 - c1 * p1) / sqrt(1 - pow(c1,2)) / sqrt(1 - pow(p1,2));
}

double g5(double p1, double p6, double c5) {
  return finv( (f(p6) * sqrt(1-pow(c5,2)) - c5 * p1) / sqrt(1-pow(p1,2)));
}

NumericVector g3(double p1, NumericVector p4, NumericVector c2,
                 NumericVector c3, NumericVector c4) {
  return (p4 * sqrt(1-pow(c2,2)) * sqrt(1 - pow(p1,2) * (1-pow(c3,2))) + c2 * sqrt(1-pow(c3,2)) * pow(p1,2)) / sqrt(1-pow(c4,2));
}


NumericVector gf8(double p1, double p2, NumericVector c7) {
  return (c7 / sqrt(1 - pow(c7,2)) * sqrt(1 - pow(p1,2)) + c7 * p1 * p2) / sqrt(1-pow(p2,2)) / sqrt(1 - pow(p1,2) * (1-pow(c7,2)));
}

double g7(double p2, double p5, double c6) {
  return finv((f(c6) * sqrt(1-pow(p5,2)) - p2 * p5) / sqrt(1-pow(p2,2)));
}


// helper method to keep main method cleaner
bool check_p2(double p1, double p2, double c6,
              NumericVector c7, NumericVector b7,
              double lb_p5, double ub_p5,
              double lb_p7s, double ub_p7s,
              int num_z) {
  bool found = false;
  
  double lb_p7c;
  double ub_p7c;
  
  if (is_false(any(is_na(c7)))) {
    NumericVector f_p8 = gf8(p1, p2, c7);
    ub_p7c = sqrt(max(b7 * pow(f_p8,2) / (1+pow(f_p8,2))));
    lb_p7c = -ub_p7c;
  } else {
    lb_p7c = -1;
    ub_p7c = 1;
  }
  
  double lb_p7 = std::max(lb_p7c, lb_p7s);
  double ub_p7 = std::min(ub_p7c, ub_p7s);
  
  if (lb_p7 <= ub_p7) {
    for (int m = 0; m < num_z && !found; m++) {
      double p5 = lb_p5 + m * (ub_p5 - lb_p5) / (num_z - 1);
      double p7 = g7(p2, p5, c6);
      found = lb_p7 <= p7 && p7 <= ub_p7;
    }
  }
  
  return found;
}



// [[Rcpp::export]]
List grid_search_cpp(int num_x, int num_y, int num_z, bool full_grid,
                     double lb_p1s, double ub_p1s,
                     double lb_p2s, double ub_p2s,
                     double lb_p3s, double ub_p3s,
                     double lb_p6s, double ub_p6s,
                     double lb_p7s, double ub_p7s,
                     NumericVector lb_p4, NumericVector ub_p4,
                     bool exist_comp_uy_bound,
                     bool exist_tsls_bound,
                     double c1,
                     NumericVector c2,
                     NumericVector c3,
                     NumericVector c4,
                     double c5,
                     double c6,
                     NumericVector c7,
                     NumericVector b7) {
  
  //Rcpp::Rcout << "Hello from C++!" << std::endl;
  
  int b_mat_dim2 = (full_grid) ? num_y : 2;
  NumericVector p1_seq(num_x);
  NumericMatrix p2_mat(num_x, b_mat_dim2);
  
  for (int i = 0; i < num_x; i++) {
    p1_seq[i] = lb_p1s + i * (ub_p1s - lb_p1s) / (num_x - 1);
    
    double lb_p2 = lb_p2s;
    double ub_p2 = ub_p2s;
    
    // pushing forward comparative OLS bounds
    if (exist_comp_uy_bound) {

      double lb_p3c;
      double ub_p3c;
      
      if (is_false(any(is_na(c2)))) {
        lb_p3c = max(g3(p1_seq[i], lb_p4, c2, c3, c4));
        ub_p3c = min(g3(p1_seq[i], ub_p4, c2, c3, c4));
      } else {
        lb_p3c = -1;
        ub_p3c = 1;
      }
      
      double lb_p3 = std::max(lb_p3c, lb_p3s);
      double ub_p3 = std::min(ub_p3c, ub_p3s);
      
      lb_p2 = std::max(lb_p2, g2(p1_seq[i], lb_p3, c1));
      ub_p2 = std::min(ub_p2, g2(p1_seq[i], ub_p3, c1));
    }
    
    
    if (lb_p2 >= ub_p2) {
      NumericMatrix::Row row = p2_mat.row(i);
      std::fill(row.begin(), row.end(), NA_REAL);
    } else if ((lb_p2 < ub_p2) && !exist_tsls_bound) {
      // b_mat_dim2 = 2 if !full_grid
      for (int j = 0; j < b_mat_dim2; j++) {
        p2_mat(i, j) = lb_p2 + j * (ub_p2 - lb_p2) / (b_mat_dim2 - 1);
      }
    } else { // There are TSLS constraints
      double lb_p5 = g5(p1_seq[i], lb_p6s, c5);
      double ub_p5 = g5(p1_seq[i], ub_p6s, c5);
      
      // computational efficiency for !full_grid case
      if (full_grid) {
        for (int k = 0; k < num_y; k++) {
          double p2 = lb_p2 + k * (ub_p2 - lb_p2) / (num_y - 1);
          bool p2_valid = check_p2(p1_seq[i], p2, c6, c7, b7,
                                   lb_p5, ub_p5, lb_p7s, ub_p7s, num_z);
          p2_mat(i, k) = p2_valid ? p2 : NA_REAL;
        }
      } else {
        bool found_lower = false;
        for (int k = 0; k < num_y && !found_lower; k++) {
          double p2 = lb_p2 + k * (ub_p2 - lb_p2) / (num_y - 1);
          found_lower = check_p2(p1_seq[i], p2, c6, c7, b7,
                                 lb_p5, ub_p5, lb_p7s, ub_p7s, num_z);
          if (found_lower) {
            p2_mat(i, 0) = p2;
          }
        }
        
        bool found_upper = false;
        for (int k = num_y - 1; k >= 0 && !found_upper; k--) {
          double p2 = lb_p2 + k * (ub_p2 - lb_p2) / (num_y - 1);
          found_upper = check_p2(p1_seq[i], p2, c6, c7, b7,
                                 lb_p5, ub_p5, lb_p7s, ub_p7s, num_z);
          if (found_upper) {
            p2_mat(i, 1) = p2;
          }
        }
        
        if (!found_lower || !found_upper) {
          p2_mat(i, 0) = NA_REAL;
          p2_mat(i, 1) = NA_REAL;
        }
        
      }
      
    }
  }
  
  return List::create(_["p1_seq"] = p1_seq, _["p2_mat"] = p2_mat);
}