#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <iostream>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
void hello_world() {
  Rcpp::Rcout << "Hello World!" << std::endl;  
}

// [[Rcpp::export]]
int get_len(vec x) {
    return x.n_elem;
}

// [[Rcpp::export]]
double compute_pdf(vec x, vec mean, mat cov) {
    //return the pdf of a multivariate Gaussian
    // test: pass
    int d = x.n_elem;
    double pdf = pow(2 * M_PI, -d * 0.5) * pow(det(cov), -0.5) * exp(- 0.5 * ((x-mean).t() * inv(cov) * (x-mean)).eval()(0, 0));

    return pdf;
}

// [[Rcpp::export]]
double compute_mixture_pdf(vec x, Rcpp::List params, vec weights) {
    // check if params if of the correct class
    // tested: pass

    double out = 0;
    Rcpp::List means = params["means"];
    Rcpp::List covs = params["covs"];

    for (int i = 0; i < weights.n_elem; i++) {
        out += weights(i) * compute_pdf(x, means[i], covs[i]);
    }

    return out;
}


// [[Rcpp::export]]
vec compute_df(vec x, vec mean, mat cov) {
    // return the gradient of a multivariate gaussian
    return -compute_pdf(x, mean, cov) * inv(cov) * (x - mean);
}


// [[Rcpp::export]]
vec compute_dUdq_mixture(vec x, Rcpp::List params, vec weights) {
    // return the gradient vector of U w.r.t x
    // target is a mixture of gaussian
    Rcpp::List means = params["means"];
    Rcpp::List covs = params["covs"];

    double pdf = compute_mixture_pdf(x, params, weights);
    vec g(size(x), fill::value(0.0));
    for (int i = 0; i < weights.n_elem; i++) {
        g -= weights[i] * compute_df(x, means[i], covs[i]);
    }
    
    return g / pdf;
}


// [[Rcpp::export]]
std::list<std::vector<double>> get_path_test() {
    std::list<std::vector<double>> vec;
    vec.push_back({1, 2});

    return vec;
}


// [[Rcpp::export]]
std::list<arma::vec> get_path(vec x, vec p, Rcpp::List params, vec weights) {

    std::list<arma::vec> path;
    int num_steps = params["num_steps"];
    double step_size = params["step_size"];
    double alpha = params["alpha"];
    double temp = alpha;
    vec point(3);


    for (int i = 0; i < num_steps; i++) {
        if (i < num_steps/2) {
            p *= sqrt(alpha);
            temp *= sqrt(alpha);
            p -= 0.5 * step_size * compute_dUdq_mixture(x, params, weights);
            x += step_size * p;
            p *= sqrt(alpha);
            temp *= sqrt(alpha);
        }
        else {
            p /= sqrt(alpha);
            temp /= sqrt(alpha);
            p -= 0.5 * step_size * compute_dUdq_mixture(x, params, weights);
            x += step_size * p;
            p -= 0.5 * step_size * compute_dUdq_mixture(x, params, weights);
            p /= sqrt(alpha);
            temp /= sqrt(alpha);
        }

        for (int j = 0; j < 2; j++) {
            point[j] = x[j];
        }
        point[2] = temp;

        path.push_back(point);
    }

    return path;
}

// [[Rcpp::export]]
double compute_acceptance(vec x0, vec p0, vec x1, vec p1, Rcpp::List params, vec weights) {
    double pdf0 = compute_mixture_pdf(x0, params, weights);
    double pdf1 = compute_mixture_pdf(x1, params, weights);
    mat mass = params["mass"];
    mat inv_mass = inv(mass);
    double dU, dK, dH;

    if (pdf0 == 0) {
        return 1;
    }
    else if (pdf1 == 0) return 0;
    else {
        dU = log(pdf1) - log(pdf0);
        dK = -0.5 * (p1.t() * inv_mass * p1).eval()(0, 0) + 0.5 * (p0.t() * inv_mass * p0).eval()(0, 0);
        dH = dU + dK;
    }

    if (dH > 0) return 1;
    else return exp(dH);
}

// [[Rcpp::export]]
vec get_next_sample(vec x, Rcpp::List params, vec weights) {
    
    int num_steps = params["num_steps"];
    double step_size = params["step_size"];
    double alpha = params["alpha"];
    mat mass = params["mass"];
    vec mean(size(x), fill::zeros);
    vec p = mvnrnd(mean, mass);
    vec p1 = vec(p);
    vec x1 = vec(x);

    for (int i = 0; i < num_steps; i++) {
        if (compute_mixture_pdf(x1, params, weights) == 0) {
            return x;
        }
        else if (i < num_steps/2) {
            p1 *= sqrt(alpha);
            p1 -= 0.5 * step_size * compute_dUdq_mixture(x1, params, weights);
            x1 += step_size * p1;
            p1 -= 0.5 * step_size * compute_dUdq_mixture(x1, params, weights);
            p1 *= sqrt(alpha);
        }
        else {
            p1 /= sqrt(alpha);
            p1 -= 0.5 * step_size * compute_dUdq_mixture(x1, params, weights);
            x1 += step_size * p1;
            p1 -= 0.5 * step_size * compute_dUdq_mixture(x1, params, weights);
            p1 /= sqrt(alpha);
        }
    }

    std::random_device rd;
    std::uniform_real_distribution<double> distr(0.0, 1.0); 
    if (compute_acceptance(x, p, x1, p1, params, weights) < distr(rd)) {
        return x;
    }
    else {
        return x1;
    }
}


// [[Rcpp::export]]
std::list<arma::vec> get_samples(vec x0, Rcpp::List params, vec weights, int n_samples) {
    std::list<arma::vec> samples;
    samples.push_back(x0);
    for (int i = 1; i < n_samples; i++) {
        samples.push_back(get_next_sample(samples.back(), params, weights));
    }

    return samples;
}

// // [[Rcpp::export]]
// std::list<arma::vec>> get_path_improved(vec x, vec p, Rcpp::List params, vec weights) {


// [[Rcpp::export]]
double compute_log_pdf(vec x, vec mean, mat cov) {
    int d = x.n_elem;
    return - (d / 2.0) * log(2 * M_PI) - log(det(cov)) - 0.5 * ((x - mean).t() * inv(cov) * (x - mean)).eval()(0, 0);
}


// After compile, this function will be immediately called using
// the below snippet and results will be sent to the R console

/*** R
hello_world() 
*/