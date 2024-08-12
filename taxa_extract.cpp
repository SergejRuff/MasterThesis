#include <Rcpp.h>
#include <boost/regex.hpp>
#include <string>
#include <vector>
#include <thread>
#include <algorithm>

// Helper function to process a chunk of data
void process_chunk(const Rcpp::CharacterVector& input, const boost::regex& re, Rcpp::CharacterVector& result, R_xlen_t start, R_xlen_t end) {
  for (R_xlen_t i = start; i < end; ++i) {
    const std::string& str = Rcpp::as<std::string>(input[i]);
    
    boost::smatch match;
    if (boost::regex_search(str, match, re) && !match.empty()) {
      result[i] = match.str(0); // Extract the first match
    } else {
      result[i] = NA_STRING; // No match found
    }
  }
}

// [[Rcpp::export]]
Rcpp::CharacterVector taxa_extract(const Rcpp::CharacterVector& input, const std::string& pattern, int num_threads = 1) {
  // Create a boost regex object
  boost::regex re(pattern);
  
  // Prepare the result vector
  Rcpp::CharacterVector result(input.size());
  
  // Determine the number of threads to use
  if (num_threads <= 0) {
    num_threads = std::thread::hardware_concurrency(); // Use hardware concurrency if num_threads is not positive
  }
  num_threads = std::max(1, num_threads); // Ensure at least one thread
  
  const R_xlen_t chunk_size = (input.size() + num_threads - 1) / num_threads;
  
  // Vector to hold thread objects
  std::vector<std::thread> threads;
  
  // Launch threads to process chunks of the input
  for (int t = 0; t < num_threads; ++t) {
    R_xlen_t start = t * chunk_size;
    R_xlen_t end = std::min(start + chunk_size, input.size());
    
    threads.emplace_back(process_chunk, std::cref(input), std::cref(re), std::ref(result), start, end);
  }
  
  // Join all threads
  for (auto& th : threads) {
    if (th.joinable()) {
      th.join();
    }
  }
  
  return result;
}
