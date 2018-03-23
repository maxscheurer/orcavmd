
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
        if (element.compare("S") == 0) {
            shell_t* shell = (shell_t*)calloc(2, sizeof(shell_t));
        *nshells = 2;
                std::vector<std::vector<double>> coeff;
                coeff.push_back({2.017722288913,-0.215082324869});
                
                coeff.push_back({0.247657587627,0.215496924081});
                
                coeff.push_back({0.116937414507,0.032235784219});
                
                coeff.push_back({0.409899773543,0.243971009984});
                
                coeff.push_back({7.411313465509,-0.185102130381});
                
                coeff.push_back({0.157964175725,0.077363489765});
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
        if (element.compare("Cl") == 0) {
            shell_t* shell = (shell_t*)calloc(2, sizeof(shell_t));
        *nshells = 2;
                std::vector<std::vector<double>> coeff;
                coeff.push_back({2.846387837613,-0.278406412992});
                
                coeff.push_back({0.349368963800,0.278943077637});
                
                coeff.push_back({0.164962857497,0.041726576369});
                
                coeff.push_back({0.689103451277,0.467032020318});
                
                coeff.push_back({12.459537714390,-0.354339730457});
                
                coeff.push_back({0.265561646276,0.148096394428});
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
        if (element.compare("F") == 0) {
            shell_t* shell = (shell_t*)calloc(2, sizeof(shell_t));
        *nshells = 2;
                std::vector<std::vector<double>> coeff;
                coeff.push_back({57.233721018966,-0.888994291217});
                
                coeff.push_back({3.475503175228,1.081300920787});
                
                coeff.push_back({1.334296115475,0.405401935353});
                
                coeff.push_back({5.704760783740,2.040840943755});
                
                coeff.push_back({1.464108504213,1.299733170342});
                
                coeff.push_back({0.497085965423,0.251251069450});
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
        if (element.compare("Br") == 0) {
            shell_t* shell = (shell_t*)calloc(2, sizeof(shell_t));
        *nshells = 2;
                std::vector<std::vector<double>> coeff;
                coeff.push_back({6.487663894806,-0.970281517671});
                
                coeff.push_back({1.272445443876,0.902317838942});
                
                coeff.push_back({0.627985792375,0.063181618513});
                
                coeff.push_back({0.336339114158,0.241118311790});
                
                coeff.push_back({2.199807320655,-0.234769827303});
                
                coeff.push_back({0.165373563848,0.059116244913});
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
        if (element.compare("P") == 0) {
            shell_t* shell = (shell_t*)calloc(2, sizeof(shell_t));
        *nshells = 2;
                std::vector<std::vector<double>> coeff;
                coeff.push_back({2.296400178172,-0.236997823436});
                
                coeff.push_back({0.281862836861,0.237454667627});
                
                coeff.push_back({0.133088235673,0.035520402252});
                
                coeff.push_back({0.337223537572,0.191156238495});
                
                coeff.push_back({6.097269396589,-0.145031276394});
                
                coeff.push_back({0.129956739637,0.060615864570});
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
