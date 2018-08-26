#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>

long int secs_ ()
{
	long int	seconds;

	time ( &seconds );
	return seconds;
}
