// Average PSF metrics
struct metrics {
    double q11 = 0, q12 = 0, q22 = 0;
    double e1 = 0, e2 = 0, r2 = 0;

    metrics() = default;

    metrics(double) {} // to enable integrate()

    explicit metrics(double t11, double t12, double t22) :
        q11(t11), q12(t12), q22(t22) {
        get_ellipticities();
    }

    void get_ellipticities() {
        r2 = q11+q22;
        e1 = (q11-q22)/r2;
        e2 = 2*q12/r2;
    }

    metrics& operator *= (double norm) {
        q11 *= norm; q12 *= norm; q22 *= norm;
        e1 *= norm; e2 *= norm; r2 *= norm;
        return *this;
    }

    metrics& operator *= (const metrics& m) {
        q11 *= m.q11; q12 *= m.q12; q22 *= m.q22;
        e1 *= m.e1; e2 *= m.e2; r2 *= m.r2;
        return *this;
    }

    metrics& operator /= (double norm) {
        q11 /= norm; q12 /= norm; q22 /= norm;
        e1 /= norm; e2 /= norm; r2 /= norm;
        return *this;
    }

    metrics& operator += (const metrics& m) {
        q11 += m.q11; q12 += m.q12; q22 += m.q22;
        e1 += m.e1; e2 += m.e2; r2 += m.r2;
        return *this;
    }

    metrics& operator -= (const metrics& m) {
        q11 -= m.q11; q12 -= m.q12; q22 -= m.q22;
        e1 -= m.e1; e2 -= m.e2; r2 -= m.r2;
        return *this;
    }

    void reset() {
        operator*=(0.0);
    }
};

metrics operator* (metrics m1, const metrics& m2) {
    return m1 *= m2;
}
metrics operator* (metrics m, double norm) {
    return m *= norm;
}
metrics operator* (double norm, metrics m) {
    return m *= norm;
}
metrics operator+ (metrics m1, const metrics& m2) {
    return m1 += m2;
}
metrics operator- (metrics m1, const metrics& m2) {
    return m1 -= m2;
}

metrics sqrt(metrics m) {
    m.q11 = sqrt(m.q11);
    m.q12 = sqrt(m.q12);
    m.q22 = sqrt(m.q22);
    m.e1 = sqrt(m.e1);
    m.e2 = sqrt(m.e2);
    m.r2 = sqrt(m.r2);
    return m;
}

namespace std {
    template<>
    struct is_arithmetic<metrics> : std::true_type {}; // to enable integrate()
}

auto get_q11 = vectorize_lambda([](const metrics& n) {
    return n.q11;
});
auto get_q12 = vectorize_lambda([](const metrics& n) {
    return n.q12;
});
auto get_q22 = vectorize_lambda([](const metrics& n) {
    return n.q22;
});
auto get_e1 = vectorize_lambda([](const metrics& n) {
    return n.e1;
});
auto get_e2 = vectorize_lambda([](const metrics& n) {
    return n.e2;
});
auto get_r2 = vectorize_lambda([](const metrics& n) {
    return n.r2;
});

// Set of average PSFs for multiple categories of galaxies
struct metrics_set {
    metrics all, qu, sf;

    void reset() {
        all.reset();
        qu.reset();
        sf.reset();
    }

    void add(uint_t it, const metrics& m) {
        all += m;
        if (it == 0) {
            qu += m;
        } else {
            sf += m;
        }
    }

    void normalize(double ntot, double nqu, double nsf) {
        all /= ntot;
        qu /= nqu;
        sf /= nsf;
    }
};

auto get_all = vectorize_lambda([](const metrics_set& n) {
    return n.all;
});
auto get_qu = vectorize_lambda([](const metrics_set& n) {
    return n.qu;
});
auto get_sf = vectorize_lambda([](const metrics_set& n) {
    return n.sf;
});

metrics_set integrate(const vec1d& z, const vec<1,metrics_set>& v,
    const vec1d& dndz, const vec1d& dndz_qu, const vec1d& dndz_sf) {

    metrics_set m;
    m.all = integrate(z, get_all(v)*dndz);
    m.qu  = integrate(z, get_qu(v)*dndz_qu);
    m.sf  = integrate(z, get_sf(v)*dndz_sf);

    return m;
}

void to_fits(std::string filename, std::string suffix, const metrics& m, const vec<1,metrics>& zm) {
    fits::update_table(filename,
        "tot_q11"+suffix, m.q11, "tot_q12"+suffix, m.q12, "tot_q22"+suffix, m.q22,
        "tot_e1"+suffix, m.e1, "tot_e2"+suffix, m.e2, "tot_r2"+suffix, m.r2
    );
    fits::update_table(filename,
        "z_q11"+suffix, get_q11(zm), "z_q12"+suffix, get_q12(zm), "z_q22"+suffix, get_q22(zm),
        "z_e1"+suffix, get_e1(zm), "z_e2"+suffix, get_e2(zm), "z_r2"+suffix, get_r2(zm)
    );
}

void to_fits(std::string filename, const metrics_set& m,
    const vec1d& z, const vec1d& dndz, const vec1d& dndz_qu, const vec1d& dndz_sf,
    const vec<1,metrics_set>& zm) {

    fits::write_table(filename, "z", z, "dndz", dndz, "dndz_qu", dndz_qu, "dndz_sf", dndz_sf);

    to_fits(filename, "",    m.all, get_all(zm));
    to_fits(filename, "_qu", m.qu,  get_qu(zm));
    to_fits(filename, "_sf", m.sf,  get_sf(zm));
}
