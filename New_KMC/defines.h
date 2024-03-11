#ifndef _DEFINES_H
#define _DEFINES_H

#define Xup 0
#define Tup 1
#define Tdn 2
#define Xdn 3
#define Ndirs 4
#define OPPdir(dir) (3-(dir))	// Opposite direction

#define COLDLAT 0
#define HOTLAT  1
#define PLUS   1
#define MINUS -1

#define FRESH  0
#define FORGET 1
#define FRESHHOT 2
#define MAXFILENAME 256

#define FORALLSITES(i,s) for(i=0,s=lattice;i<Volume;i++,s++)

#endif // _DEFINES_H
