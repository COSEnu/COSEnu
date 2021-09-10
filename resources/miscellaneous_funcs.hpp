/*---------------------------------------------------------------------------*/
template <typename T>
void string_to_type(std::string value, T &var)
{
    std::stringstream svalue(value);
    svalue >> var;
}
/*---------------------------------------------------------------------------*/
template <typename T>
void cssl_to_vec(std::string cssl, std::vector<T> &vec) //cssl -> comma-seperated-string-list
{
    int start = cssl.find_first_of("[");
    int end = cssl.find_first_of("]");
    cssl = cssl.substr(start + 1, end - 1);

    // In case any accidantal space in the begining of cssl.
    cssl = cssl.substr(cssl.find_first_not_of(" "), end - 1) + ",;";
    std::string c = "";

    while (true)
    {
        c = cssl.substr(0, cssl.find_first_of(","));
        cssl = cssl.substr(cssl.find_first_of(","), cssl.find_first_of("]"));
        cssl = cssl.substr(cssl.find_first_not_of(" ") + 1, cssl.find_first_of("]"));
        std::stringstream svalue(c);
        T tval;
        svalue >> tval;
        vec.push_back(tval);

        if (cssl == ";")
            break;
    }
}
/*---------------------------------------------------------------------------*/
std::string draw(const int l, const string c)
{
    // write a string l times on to std::out
    std::string s = "";
    for (int i = 0; i < l; i++)
    {
        s += c;
    }
    return s;
}
/*---------------------------------------------------------------------------*/

template <typename T>
int sign(T val)
{
    return (T(0) < val) - (val < T(0));
}
/*---------------------------------------------------------------------------*/

inline void swap(FieldVar **a, FieldVar **b)
{
    FieldVar *tmp = *a;
    *a = *b;
    *b = tmp;
}

/*---------------------------------------------------------------------------*/

double g(double v, double v0, double sigma)
{
    double exponant = (v - v0) * (v - v0) / (2.0 * sigma * sigma);
    double N = sigma * sqrt(M_PI / 2.0) * (erf((1.0 + v0) / sigma / sqrt(2.0)) + erf((1.0 - v0) / sigma / sqrt(2.0)));
    return exp(-exponant) / N;
}

/*---------------------------------------------------------------------------*/

inline double eps(double z, double z0)
{
    return 1e-1 * exp(-(z - z0) * (z - z0) / 50.0);
}

/*---------------------------------------------------------------------------*/

inline double eps_(double z, double z0)
{
    double e = eps(z, z0);
    return sqrt(1.0 - e * e);
}

/*---------------------------------------------------------------------------*/

inline double random_amp(double a)
{
    return a * rand() / RAND_MAX;
}
/*---------------------------------------------------------------------------*/

double hev(double x, double x0)
{
    /* Heviside theta function */
    return 0.5 * (1.0 + sign(x - x0));
}

/*---------------------------------------------------------------------------*/

double gauss(double x, double x0, double sigma)
{
    double exponant = (x - x0) * (x - x0) / (2.0 * sigma * sigma);
    return exp(-exponant);
}

/*---------------------------------------------------------------------------*/

double L(double x, unsigned int order)
{
    // Legendre polynomials.
    double value;
    switch (order)
    {
    case 0:
        value = 1;
        break;
    case 1:
        value = x;
        break;
    case 2:
        value = 0.5*(3*x*x - 1);
        break;
    case 3:
        value = 0.5*(5*pow(x, 3)-3*x);
        break;
    default:
        value = 0;
    }
    return (value);
}

/*---------------------------------------------------------------------------*/
