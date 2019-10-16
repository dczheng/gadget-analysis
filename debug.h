//#define TURNOFF_ALL_TEST

#ifndef TURNOFF_ALL_TEST

//#define GROUP_POT_TEST
//#define OUTPUTGROUP_DEBUG

#ifdef GROUP_POT_TEST2
#ifndef GROUP_POT_TEST
#define GROUP_POT_TEST
#endif
#endif

//#define FOF_TEST
//#define MATH_TEST
//#define QTRAP_DEBUB
//#define PART_RADIO_TEST
#define PART_RADIO_TEST2
//#define PART_RADIO_TEST3

#ifdef PART_RADIO_TEST
#ifndef RAD
#define RAD
#endif
#ifndef RADLARMOR
#define RADLARMOR
#endif
#endif

#endif
