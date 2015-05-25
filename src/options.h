#pragma once

#include <map>

#include <getopt.h>

using namespace std;

namespace options
{
    option long_options[] = {
            {"homo-only", no_argument, nullptr, 0},
            {"hetero-only", no_argument, nullptr, 0},
            {"binding-cutoff", required_argument, nullptr, 0 },
            {"nonbinding-cutoff", required_argument, nullptr, 0}
    };
    map<string, string> opt_map;
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

