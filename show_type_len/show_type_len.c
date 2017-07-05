#include "stdio.h"
void main() {
    char *buf = "%20s : %i\n";
    printf( buf, "char", sizeof( char ) );
    printf( buf, "unsigned char", sizeof( unsigned char ) );
    printf( buf, "short", sizeof( short ) );
    printf( buf, "unsigned short", sizeof( unsigned short ) );
    printf( buf, "int", sizeof( int ) );
    printf( buf, "unsigend int", sizeof( signed int ) );
    printf( buf, "long", sizeof( long ) );
    printf( buf, "unsigned long", sizeof( unsigned long ) );
    printf( buf, "long long", sizeof( long long ) );
    printf( buf, "unsigned long long", sizeof( unsigned long long ) );
    printf( buf, "float", sizeof( float ) );
    printf( buf, "double", sizeof( double ) );
    printf( buf, "long double", sizeof( long double ) );
}
