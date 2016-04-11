#include "interfaces.hpp"
#include "CC_likelihood.hpp"

int main()
{
    Settings settings;

    settings.nDims         = 20;
    settings.nDerived      = 1;
    settings.nlive         = 500;
    settings.num_repeats   = settings.nDims*5;
    settings.do_clustering = false;

    settings.precision_criterion = 1e-3;

    settings.base_dir      = "chains";
    settings.file_root     = "test";

    settings.write_resume  = false;
    settings.read_resume   = false;
    settings.write_live    = true;
    settings.write_dead    = false;
    settings.write_stats   = false;

    settings.equals        = false;
    settings.posteriors    = false;
    settings.cluster_posteriors = false;

    settings.feedback      = 1;
    settings.update_files  = settings.nlive;

    settings.boost_posterior= 5.0;

    setup_loglikelihood();
    run_polychord(loglikelihood, prior, settings) ;


}
