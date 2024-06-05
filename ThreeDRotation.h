//modified from http://www.programming-techniques.com/2012/03/3d-rotation-algorithm-about-arbitrary.html
//there is typos in the matrix formula in the blog above
// the correct one can be found here https://docs.google.com/viewer?a=v&pid=sites&srcid=ZGVmYXVsdGRvbWFpbnxnbGVubm11cnJheXxneDoyMTJiZTZlNzVlMjFiZTFi



#ifndef _3DROTATE_
#define _3DROTATE_

void PerformRotation(double point_x, double point_y, double point_z,//point to rotate
	                        double u_axis, double v_axis, double w_axis,//axis to rotate about
	                        double a, double b, double c, //point the axis of rotation passes through
	                        double angle,//angle of rotation
	                        double&rot_point_x, double&rot_point_y, double&rot_point_z);//point after roation

#endif /*_3DROTATE_*/


