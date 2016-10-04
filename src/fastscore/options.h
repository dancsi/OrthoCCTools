#pragma once

#include <map>
#include <cstdio>
#include <sstream>

#ifdef _MSC_VER
#include "getopt.h"
#else
#include <getopt.h>
#endif

using namespace std;

namespace options
{

    option long_options[] = {
            {"out-name", required_argument, nullptr, 0},
            {"align", required_argument, nullptr, 0}
    };
    map<string, string> opt_map;
    void usage(char** argv)
    {
        fprintf(stderr, "USAGE: %s INPUT [OPTIONS...]\n", argv[0]);
        fprintf(stderr, "Available options:\n");
        for(int i=0;i<sizeof(long_options)/sizeof(option);i++)
        {
            fprintf(stderr, "\t--%s%s\n", long_options[i].name, long_options[i].has_arg==required_argument?"=PARAM":"");
        }
    }
    void parse(int argc, char **argv)
    {
        int c, index;
        while((c = getopt_long(argc, argv, "", long_options, &index))!=-1)
        {
            string name = string(long_options[index].name);
            switch(c)
            {
                case 0:
                    if(long_options[index].has_arg == no_argument) opt_map.insert(make_pair(name, string("1")));
                    else if(long_options[index].has_arg == required_argument) opt_map.insert(make_pair(name, string(optarg)));
                    break;
                default:
                    break;
            }
        }
    }

    template<typename T>
    T get(string name, T default_ret)
    {
        auto it = opt_map.find(name);
        T ret;
        if(it!=opt_map.end())
        {
            stringstream ss(it->second);
            ss>>ret;
        }
        else
        {
            ret = default_ret;
        }
        cerr<<"Value of "<<name<<" is "<<ret<<endl;
        return ret;
    }

};

