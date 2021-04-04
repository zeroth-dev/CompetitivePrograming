
#include <iostream>
#include <string>
#include <sstream>
#include <curlpp/cURLpp.hpp>
#include <curlpp/Options.hpp>
#include <curlpp/Easy.hpp>
#include <curlpp/Infos.hpp>

// RAII cleanup
std::string uniprotDownload(const char* id){
// Send request and get a result.
// Here I use a shortcut to get it in a string stream ...

  //  curlpp::Cleanup cleaner;
std::stringstream response;
curlpp::options::Url myUrl(id); // Don't cast to string
curlpp::Easy myRequest;
myRequest.setOpt(curlpp::options::WriteStream(&std::cout));
myRequest.setOpt(myUrl);
myRequest.setOpt( new curlpp::options::WriteStream( &response ) );
myRequest.setOpt( new curlpp::options::FollowLocation(1));

myRequest.perform();

// There is some shorcut within curlpp that allow you to write shorter code
// like this:
return response.str();
}