//
// Created by Yuzhen Feng on 19/2/2025.
//

#ifndef COST_H
#define COST_H

#include "utils.h"

struct bpr {
    double capacity;
    double fft;
    double alpha;
    double beta;

    double betap1;
    double betam1;
    double constant;

    bpr() : capacity(0.), fft(0.), alpha(0.), beta(0.), betap1(0.), betam1(0.), constant(0.) {}

    double operator()(const double& flow) const {
        return fft * (1. + alpha * math2::opt_pow(flow / capacity, beta));
    }

    [[nodiscard]] double derivative(const double& flow) const {
        return beta * constant * math2::opt_pow(flow / capacity, betam1);
    }

    [[nodiscard]] double integral(const double& flow) const {
        return flow * fft + (fft * alpha * capacity / betap1) * math2::opt_pow(flow / capacity, betap1);
    }

    void initialize(const double& _capacity, const double& _fft=0, const double& _alpha=0.15, const double& _beta=4) {
        this->capacity = _capacity;
        this->fft = _fft;
        this->alpha = _alpha;
        this->beta = _beta;

        this->betap1 = this->beta + 1;
        this->betam1 = this->beta - 1;
        this->constant = (fft * alpha) / capacity;
    }

    void update(const double& flow, double& cost, double& derivative) const {
        cost = this->operator()(flow);
        derivative = this->derivative(flow);
    }
};

struct linear {
    double vot;
    double capacity;

    linear() : vot(0.), capacity(0.) {}

    void initialize(const double& _capacity, const double& _vot=1) {
        this->vot = _vot;
        this->capacity = _capacity;
    }

    double operator()(const double& flow) const {
        return flow / capacity;
    }

    [[nodiscard]] double derivative(const double& flow) const {
        return 1. / capacity;
    }

    [[nodiscard]] double integral(const double& flow) const {
        return 0.5 * flow * flow / capacity;
    }

    void update(const double& flow, double& cost, double& derivative) const {
        cost = this->operator()(flow);
        derivative = this->derivative(flow);
    }
};

struct quad {
    double c;
    double d;

    quad() : c(0.), d(0.) {}

    void initialize(const double& _c, const double& _d) {
        this->c = _c;
        this->d = _d;
    }

    double operator()(const double& flow) const {
        return 2 * c * flow + d;
    }

    [[nodiscard]] double derivative(const double& flow) const {
        return 2 * c;
    }

    [[nodiscard]] double integral(const double& flow) const {
        return c * flow * flow + d * flow;
    }

    void update(const double& flow, double& cost, double& derivative) const {
        cost = this->operator()(flow);
        derivative = this->derivative(flow);
    }
};

#endif //COST_H
