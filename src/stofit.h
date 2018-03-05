
#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <cctype>
#include <sstream>
#include "qmplugin.h"

shell_t* get_bas_4element(std::string element, int* nshells) {
        prim_t* prim;


        if (element.compare("O") == 0) {
            shell_t* shell = (shell_t*)calloc(2, sizeof(shell_t));
        *nshells = 2;
                std::vector<std::vector<double>> coeff;
                coeff.push_back({37.209476655259,-0.643650567889});

                coeff.push_back({2.259536021100,0.782884613096});

                coeff.push_back({0.867468675390,0.293519529306});

                coeff.push_back({5.248151505256,1.838739966806});

                coeff.push_back({1.346921201700,1.171022824588});

                coeff.push_back({0.457299185114,0.226370107144});
                for (size_t ns = 0; ns < 2; ++ns) {
                    prim = (prim_t*) calloc(3, sizeof(prim_t));
                    shell[ns].prim = prim;
                    shell[ns].numprims = 3;
                    shell[ns].type = ns;

                    for (size_t nbas = 0; nbas < shell[ns].numprims; ++nbas) {
                        std::vector<double> block = coeff[3*ns+nbas];
                        prim[nbas].exponent = block[0];
                        prim[nbas].contraction_coeff = block[1];
                    }
            }

            return shell;        }
        if (element.compare("C") == 0) {
            shell_t* shell = (shell_t*)calloc(2, sizeof(shell_t));
        *nshells = 2;
                std::vector<std::vector<double>> coeff;
                coeff.push_back({6.323427520939,-0.170362652855});

                coeff.push_back({0.383988530469,0.207215383968});

                coeff.push_back({0.147418770394,0.077689305614});

                coeff.push_back({3.120109517380,0.959897054268});

                coeff.push_back({0.800766070944,0.611321546328});

                coeff.push_back({0.271871636773,0.118174403637});
                for (size_t ns = 0; ns < 2; ++ns) {
                    prim = (prim_t*) calloc(3, sizeof(prim_t));
                    shell[ns].prim = prim;
                    shell[ns].numprims = 3;
                    shell[ns].type = ns;

                    for (size_t nbas = 0; nbas < shell[ns].numprims; ++nbas) {
                        std::vector<double> block = coeff[3*ns+nbas];
                        prim[nbas].exponent = block[0];
                        prim[nbas].contraction_coeff = block[1];
                    }
            }

            return shell;        }
        if (element.compare("N") == 0) {
            shell_t* shell = (shell_t*)calloc(2, sizeof(shell_t));
        *nshells = 2;
                std::vector<std::vector<double>> coeff;
                coeff.push_back({10.618247577473,-0.251304127039});

                coeff.push_back({0.644790387797,0.305666061806});

                coeff.push_back({0.247544388929,0.114600487844});

                coeff.push_back({4.920990456194,1.696594342929});

                coeff.push_back({1.262956370861,1.080495739203});

                coeff.push_back({0.428791913366,0.208870340625});
                for (size_t ns = 0; ns < 2; ++ns) {
                    prim = (prim_t*) calloc(3, sizeof(prim_t));
                    shell[ns].prim = prim;
                    shell[ns].numprims = 3;
                    shell[ns].type = ns;

                    for (size_t nbas = 0; nbas < shell[ns].numprims; ++nbas) {
                        std::vector<double> block = coeff[3*ns+nbas];
                        prim[nbas].exponent = block[0];
                        prim[nbas].contraction_coeff = block[1];
                    }
            }

            return shell;        }
        if (element.compare("H") == 0) {
            shell_t* shell = (shell_t*)calloc(1, sizeof(shell_t));
        *nshells = 1;
                std::vector<std::vector<double>> coeff;
                coeff.push_back({2.086464404671,0.190952726016});

                coeff.push_back({0.380059714778,0.184680090127});

                coeff.push_back({0.102859818272,0.057556211903});
                for (size_t ns = 0; ns < 1; ++ns) {
                    prim = (prim_t*) calloc(3, sizeof(prim_t));
                    shell[ns].prim = prim;
                    shell[ns].numprims = 3;
                    shell[ns].type = ns;

                    for (size_t nbas = 0; nbas < shell[ns].numprims; ++nbas) {
                        std::vector<double> block = coeff[3*ns+nbas];
                        prim[nbas].exponent = block[0];
                        prim[nbas].contraction_coeff = block[1];
                    }
            }

            return shell;        }
 }
