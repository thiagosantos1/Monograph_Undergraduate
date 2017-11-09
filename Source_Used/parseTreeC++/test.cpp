#include <iostream>
#include <string>
#include <regex>
#include <iterator>
#include <cstring>
#include <string.h>

using namespace std;
int main ()
{

  std::string s ("(((A:0.3,B:0.4)C:0.3,D:0.4)E:0.2,(F:0.3,(G:0.3,H:0.4)I:0.4)J:0.5)K;");

  char cs[67] = {'(','(','(','A',':','0','.','3',',','B',':','0','.','4',')','C',':','0',
                  '.','3',',','D',':','0','.','4',')','E',':','0','.','2',',','(',
                  'F',':','0','.','3',',','(','G',':','0','.','3',',','H',':','0','.',
                  '4',')','I',':','0','.','4',')','J',':','0','.','5',')','K',';'};

  std::regex e ("\\([^(^)]+\\)");
  std::smatch m;
  int i=0;

  while(i<4){

    

    string test = std::regex_replace (cs,e,"");
    const char * tes;

    tes = test.c_str();
    
    memset(cs, 0, sizeof cs);
    strncpy (cs,tes,sizeof tes);

    cout<<tes<<endl;

    
    

    i++;

  }


  // using string/c-string (3) version:
  //std::cout << std::regex_replace (s,e,"");

/*
  // using range/c-string (6) version:
  std::string result;
  std::regex_replace (std::back_inserter(result), s.begin(), s.end(), e, "$2");
  std::cout << result;

  // with flags:
  std::cout << std::regex_replace (s,e,"$1 and $2",std::regex_constants::format_no_copy);
  std::cout << std::endl;
*/
  return 0;
}