#include "dataio.h"
#include <fstream>
#include <filesystem>

using namespace std;

std::vector<int> splitStringToInt( const std::string& s )
{
	vector<int> vd;
	size_t i,j;
	i = 0; j = s.find(',');
	while (i != string::npos)
	{
		if (i < j)
			vd.push_back(atoi(s.substr(i,j-i).c_str()));

		if (j == string::npos) i = j;
		else {
			i = j + 1;
			j = s.find(',',j+1);
		}
	}
	return vd;
}

std::vector<double> splitStringToDouble( const std::string& s )
{
	vector<double> vd;
	size_t i,j;
	i = 0; j = s.find(',');
    while (i != string::npos)
	{
		if (i < j)
			vd.push_back(atof(s.substr(i,j-i).c_str()));

		if (j == string::npos) i = j;
		else {
			i = j + 1;
			j = s.find(',',j+1);
		}
	}
	return vd;
}

bool fileExist(const std::string& filename)
{
    std::tr2::sys::path myfile(filename);
    return std::tr2::sys::exists(myfile);
}



std::vector<std::string> readFileList(const std::string &filename)
{
    if (!fileExist(filename)) {
        throw runtime_error("Cannot open list file!");
    }

    vector<string> result;
    ifstream list_file(filename);
    while (!list_file.eof()) {
        string line_content;
        getline(list_file, line_content);
        for (auto iter = line_content.begin(); iter != line_content.end();) {
            if (*iter == ' ' || *iter == '\t') iter = line_content.erase(iter);
            else ++iter;
        }
        if (line_content == "") continue;
        if (line_content[0] == '#') continue;
        result.push_back(line_content);
    }
    return result;
}

std::pair<std::string, std::string> splitFileName(const std::string& file_name)
{
    size_t dot_pos = file_name.rfind('.'), slashPos = file_name.rfind('/');
    std::string name = file_name.substr(slashPos + 1, dot_pos - slashPos - 1);
    std::string ext = file_name.substr(dot_pos, file_name.size() - dot_pos);
    return make_pair(name, ext);
}
