#ifndef __HllRiemanSolver_HEADER__
#define __HllRiemanSolver_HEADER__

class RiemannSolver
{

} ;


class HllRiemannSolver : public RiemannSolver
{
public:
  int IntercellFlux(const double *pl, const double *pr, double *U,
		    double *F, double s, int dim);

  int IntercellFlux(const double *pl, const double *pr,
		    const double *ul, const double *ur,
		    const double *fl, const double *fr,
		    double *U,
		    double *F,
		    double s,
		    int dim);
} ;


#endif // __HllRiemanSolver_HEADER__
