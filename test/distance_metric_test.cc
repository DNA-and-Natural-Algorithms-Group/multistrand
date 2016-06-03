#include <iostream>
#include <string.h>
#include <vector>

using namespace std;

typedef std::vector<int> intvec;


int distancemetric( char *our_struc, char *stop_struc )
{
  int loop,len;
  intvec our_pairs, stop_pairs;
  int remaining_distance = 0;

  len = strlen(our_struc);
  if( len != strlen(stop_struc) )
	return -1;  // something weird happened, as it should have the
				   // same ID list...

  //  pairqueue = new int [len+2];
  //  pairqueue.resize( len+2, -1 ); // default entry is -1, unpaired

  for( loop = 0; loop < len; loop ++ )
    {
      if( stop_struc[loop] != '*' )
        if( our_struc[loop] != stop_struc[loop] )
          {
			remaining_distance--;
          }

      if( our_struc[loop] == '(')
		our_pairs.push_back( loop );
	  if( stop_struc[loop] == '(' )
		stop_pairs.push_back( loop );

      if( our_struc[loop] == ')' && stop_struc[loop] == ')' )
        {
		  if( our_pairs.back() != stop_pairs.back() )
			{
			  remaining_distance--; // for position loop, which had
									// ),) but they were paired wrong
			  if( our_struc[stop_pairs.back()] == '(' )
				remaining_distance--; // for the position we were
									  // paired with in stop_struc,
									  // because it was ( in our_struc
									  // as well, but paired wrong
			}
		  our_pairs.pop_back();
		  stop_pairs.pop_back();
        }
      else 
		{
		  if( our_struc[loop] == ')' )
			our_pairs.pop_back();
		  if( stop_struc[loop] == ')' )
			{
			  if( our_struc[stop_pairs.back()] == '(' )
				remaining_distance--;  // for the position we were
									   // paired with in stop_struc,
									   // because it was ( in our
									   // struc but paired wrong. Note
									   // we have already subtracted
									   // for current position loop,
									   // as our_struc[loop] !=
									   // stop_struc[loop] in this
									   // conditional block.
			  stop_pairs.pop_back();
			}
		}

	  // if( remaining_distance < 0 )
	  // 	return false;
    }
  return -remaining_distance;
}


int main( int argc, char **argv )
{
  char base_struc[] =        "...(((....)))...";
  char shift_left_struc[] =  "..(.((....)))...";
  char shift_right_struc[] = "...(((....)).)..";
  char missing_struc[] =     "....((....))....";
  char added_struc[] =       "..((((....))))..";
  char shift_loose_struc[] = "**((((....)**)))";

  cout << "Base vs base:" << distancemetric(base_struc, base_struc ) << "\n";
  cout << "Base vs L shift:" << distancemetric(base_struc, shift_left_struc ) << "\n";
  cout << "Base vs R shift:" << distancemetric(base_struc, shift_right_struc ) << "\n";
  cout << "Base vs missing:" << distancemetric(base_struc, missing_struc ) << "\n";
  cout << "Base vs added:" << distancemetric(base_struc, added_struc ) << "\n";
  cout << "Base vs loose shift:" << distancemetric(base_struc, shift_loose_struc ) << "\n";


}
