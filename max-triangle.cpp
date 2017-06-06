#include "max-triangle.hpp"

#define getmod(arg1,arg2) arg1[(((arg2)+(arg1.size()))%(arg1.size()))]

bool is_convex(std::vector<mpq_class> polygon) {
    //check if polygon is a strictly convex polygon with vertices listed in CCW order

    unsigned int i;
    mpq_class x;

    for (i = 0; i<polygon.size(); i++) {
	x  = (getmod(polygon,2*i+2) - getmod(polygon,2*i+0))*(getmod(polygon,2*i+5) - getmod(polygon,2*i+1));
	x -= (getmod(polygon,2*i+4) - getmod(polygon,2*i+0))*(getmod(polygon,2*i+3) - getmod(polygon,2*i+1));
	if (x <= 0) return false;
    }

    return true;

}

struct state anchored_triangle(std::vector<mpq_class> polygon, mpq_class nx, mpq_class ny){

    struct state my_state;
    mpq_class minx,maxx,x;
    mpq_class ax,ay,bx,by,cx,cy;
    mpq_class ebx,eby,ecx,ecy;
    mpq_class dbt, dct, t0, tq;
    mpq_class ee, neb, nec;
    unsigned int i;

    if (polygon.size() < 3) { my_state.status = status_no_interior; return my_state; }
    if (!is_convex(polygon)) { my_state.status = status_not_convex; return my_state; }
    my_state.status = status_ok;

    maxx = minx = nx*polygon[0]+ny*polygon[1];
    my_state.ai = my_state.ci = 0;
    for (i = 1; i < polygon.size()/2; i++) {
	x = nx*polygon[2*i+0]+ny*polygon[2*i+1];
	if (x>maxx) { maxx=x; my_state.ci=i; }
	if (x<minx) { minx=x; my_state.ai=i; }
    }
    my_state.bi=my_state.ci-1;
    my_state.bt=1;
    my_state.ct=0;

    ax = getmod(polygon, 2*my_state.ai+0);
    ay = getmod(polygon, 2*my_state.ai+1);
    cx = getmod(polygon, 2*my_state.ci+0);
    cy = getmod(polygon, 2*my_state.ci+1);
    bx = cx;
    by = cy;
    ebx = getmod(polygon, 2*my_state.bi+2) - getmod(polygon, 2*my_state.bi+0);
    eby = getmod(polygon, 2*my_state.bi+3) - getmod(polygon, 2*my_state.bi+1);
    ecx = getmod(polygon, 2*my_state.ci+2) - getmod(polygon, 2*my_state.ci+0);
    ecy = getmod(polygon, 2*my_state.ci+3) - getmod(polygon, 2*my_state.ci+1);

    while(true) {
	neb = ebx*nx + eby*ny; 
	nec = ecx*nx + ecy*ny;
	ee = ebx*ecy - ecx*eby;
	if (neb < 0) { my_state.status = status_runtime_error; return my_state; }
	if (nec > 0) { my_state.status = status_runtime_error; return my_state; }

	if (neb == 0) { //eb.n = 0, advance b
	    my_state.bi -= 1; 
	    my_state.bt = 1;
	    bx = getmod(polygon, 2*my_state.bi+2);
	    by = getmod(polygon, 2*my_state.bi+3);
	    ebx = getmod(polygon, 2*my_state.bi+2) - getmod(polygon, 2*my_state.bi+0);
	    eby = getmod(polygon, 2*my_state.bi+3) - getmod(polygon, 2*my_state.bi+1);
	    continue;
	}
	if (nec == 0) { //ec.n = 0, advance c
	    my_state.ci += 1; 
	    my_state.ct = 0;
	    cx = getmod(polygon, 2*my_state.ci+0);
	    cy = getmod(polygon, 2*my_state.ci+1);
	    ecx = getmod(polygon, 2*my_state.ci+2) - getmod(polygon, 2*my_state.ci+0);
	    ecy = getmod(polygon, 2*my_state.ci+3) - getmod(polygon, 2*my_state.ci+1);
	    continue;
	}
	if (ee <= 0) { //eb and ec are parallel or converging, any advancement reduces area, so done
	    return my_state;
	}

	//dbt  = nx*eby*ecx*(ax-bx) + nx*ecy*ebx*(2*bx-ax-cx);
	//dbt += ny*ebx*ecy*(by-ay) + ny*ecx*eby*(ay+cy-2*by);
	//dbt /= 2*(ebx*ecy - ecx*eby)*(ebx*nx + eby*ny);

	//dct  = nx*ecy*ebx*(ax-cx) + nx*eby*ecx*(2*cx-ax-bx);
	//dct += ny*ecx*eby*(cy-ay) + ny*ebx*ecy*(ay+by-2*cy);
	//dct /= 2*(ebx*ecy - ecx*eby)*(ecx*nx + ecy*ny);
	
	tq  = nx*(by*ebx*ecx - cy*ebx*ecx - ax*eby*ecx + cx*eby*ecx + ax*ebx*ecy - bx*ebx*ecy);
	tq -= ny*(ay*eby*ecx - by*eby*ecx - ay*ebx*ecy + cy*ebx*ecy + bx*eby*ecy - cx*eby*ecy);

	dbt = -tq/(2*neb*ee);
	dct = tq/(2*nec*ee);

	//std::cout << my_state.bi << ", " << my_state.ci << ", " << my_state.bt << ", " << my_state.ct << ", " << dbt << ", " << dct << ", " << bx << ", " << by << ", " << cx << ", " << cy << std::endl;

	if (dbt <= 0 || dct <= 0) return my_state;

	if (dbt < my_state.bt && dct < 1 - my_state.ct) {
	    my_state.bt -= dbt;
	    my_state.ct += dct;
	    return my_state;
	}

	if ( (1 - my_state.ct)*dbt < my_state.bt*dct ) {
	    t0 = (1 - my_state.ct)*dbt/dct;
	    bx -= t0*ebx;
	    by -= t0*eby;
	    my_state.bt -= t0;
	    my_state.ci += 1; 
	    my_state.ct = 0;
	    cx = getmod(polygon, 2*my_state.ci+0);
	    cy = getmod(polygon, 2*my_state.ci+1);
	    ecx = getmod(polygon, 2*my_state.ci+2) - getmod(polygon, 2*my_state.ci+0);
	    ecy = getmod(polygon, 2*my_state.ci+3) - getmod(polygon, 2*my_state.ci+1);
	} else {
	    t0 = my_state.bt*dct/dbt;
	    cx += t0*ecx;
	    cy += t0*ecy;
	    my_state.ct += t0;
	    my_state.bi -= 1; 
	    my_state.bt = 1;
	    bx = getmod(polygon, 2*my_state.bi+2);
	    by = getmod(polygon, 2*my_state.bi+3);
	    ebx = getmod(polygon, 2*my_state.bi+2) - getmod(polygon, 2*my_state.bi+0);
	    eby = getmod(polygon, 2*my_state.bi+3) - getmod(polygon, 2*my_state.bi+1);
	}
    }
}


struct state maximum_triangle(std::vector<mpq_class> polygon){

    struct state my_state = anchored_triangle(polygon,1,0);
    my_state.nx = 1; my_state.ny = 0;
    std::vector<unsigned int> ret = {0,0,0};
    unsigned int sz = polygon.size()/2;
    bool max_at_vert=false;

    if (my_state.status != status_ok) {return my_state;}

    int ai_start = my_state.ai;
    mpq_class ax,ay,bx,by,cx,cy;
    mpq_class pax,pay,pbx,pby,pcx,pcy;
    mpq_class eax,eay,ebx,eby,ecx,ecy;
    mpq_class t1,t2,t3;
    mpq_class ee,fb,fc,tq;
    mpq_class area;

    ax = getmod(polygon, 2*my_state.ai+0);
    ay = getmod(polygon, 2*my_state.ai+1);
    bx = getmod(polygon, 2*my_state.bi+0);
    by = getmod(polygon, 2*my_state.bi+1);
    cx = getmod(polygon, 2*my_state.ci+0);
    cy = getmod(polygon, 2*my_state.ci+1);

    eax = getmod(polygon, 2*my_state.ai+2) - getmod(polygon, 2*my_state.ai+0);
    eay = getmod(polygon, 2*my_state.ai+3) - getmod(polygon, 2*my_state.ai+1);
    ebx = getmod(polygon, 2*my_state.bi+2) - getmod(polygon, 2*my_state.bi+0);
    eby = getmod(polygon, 2*my_state.bi+3) - getmod(polygon, 2*my_state.bi+1);
    ecx = getmod(polygon, 2*my_state.ci+2) - getmod(polygon, 2*my_state.ci+0);
    ecy = getmod(polygon, 2*my_state.ci+3) - getmod(polygon, 2*my_state.ci+1);

    //mpq_class amax = (bx-ax)*(cy-ay) - (cx-ax)*(by-ay);
    //ret[0] = my_state.ai; ret[1] = my_state.bi; ret[2] = my_state.ci;
    mpq_class amax = -1;


    while (my_state.ai <= ai_start + (int) sz) {
	if (my_state.nx*eax + my_state.ny*eay == 0) { //support line coincides with forward edge at A, advance A
	    my_state.ai += 1;
	    ax = getmod(polygon, 2*my_state.ai+0);
	    ay = getmod(polygon, 2*my_state.ai+1);
	    eax = getmod(polygon, 2*my_state.ai+2) - getmod(polygon, 2*my_state.ai+0);
	    eay = getmod(polygon, 2*my_state.ai+3) - getmod(polygon, 2*my_state.ai+1);
	    continue;
	}
	if (my_state.bt > 1) {my_state.status=status_runtime_error; return my_state;}
	if (my_state.ct > 1)  {my_state.status=status_runtime_error; return my_state;}
	if (my_state.bt == 1) {
	    my_state.bi += 1;
	    my_state.bt = 0;
	    bx = getmod(polygon, 2*my_state.bi+0);
	    by = getmod(polygon, 2*my_state.bi+1);
	    ebx = getmod(polygon, 2*my_state.bi+2) - getmod(polygon, 2*my_state.bi+0);
	    eby = getmod(polygon, 2*my_state.bi+3) - getmod(polygon, 2*my_state.bi+1);
	    continue;
	}
	if (my_state.ct == 1) {
	    my_state.ci += 1;
	    my_state.ct = 0;
	    cx = getmod(polygon, 2*my_state.ci+0);
	    cy = getmod(polygon, 2*my_state.ci+1);
	    ecx = getmod(polygon, 2*my_state.ci+2) - getmod(polygon, 2*my_state.ci+0);
	    ecy = getmod(polygon, 2*my_state.ci+3) - getmod(polygon, 2*my_state.ci+1);
	    continue;
	}
	area = (bx-ax)*(cy-ay) - (cx-ax)*(by-ay);
	if ( (area > amax) || ( (area == amax) && !max_at_vert ) ) {
	    amax = area;
	    if (my_state.ct == 0 && my_state.bt == 0) {
		max_at_vert = true;
		ret[0] = (my_state.ai+sz)%sz;
		ret[1] = (my_state.bi+sz)%sz;
		ret[2] = (my_state.ci+sz)%sz;
	    } else max_at_vert = false;
	}

	ee = ebx*ecy - ecx*eby;
	fb = ebx*(ay + by - 2*cy) - eby*(ax + bx - 2*cx);
	fc = ecx*(ay + cy - 2*by) - ecy*(ax + cx - 2*bx);
	tq  = my_state.nx*(by*ebx*ecx - cy*ebx*ecx - ax*eby*ecx + cx*eby*ecx + ax*ebx*ecy - bx*ebx*ecy);
	tq -= my_state.ny*(ay*eby*ecx - by*eby*ecx - ay*ebx*ecy + cy*ebx*ecy + bx*eby*ecy - cx*eby*ecy);
	//tq > 0 ---> move right
	//tq < 0 ---> move left
	//tq = 0 ---> swing
	
	if (tq > 0) {
	    //C stays fixed, B moves along eb
	    //possible stopping conditions: 
	    //   (1) B hits end of edge,
	    //   (2) BC becomes parallel to ea,
	    //   (3) tq becomes 0
	    t1 = 1 - my_state.bt;

	    t2 = (by*eax - cy*eax - bx*eay + cx*eay)/(eay*ebx - eax*eby);

	    t3  = ebx*ecx*(by-cy)*(by-cy);
	    t3 += eby*ecy*(bx-cx)*(bx-cx);
	    t3 += ebx*ecy*(-ay*bx + ax*by - bx*by + ay*cx - ax*cy + 2*bx*cy - cx*cy);
	    t3 += eby*ecx*( ay*bx - ax*by - bx*by - ay*cx + 2*by*cx + ax*cy - cx*cy);
	    t3 /= ee*fb;


	    if (t3 < t1) t1 = t3;
	    if (t2 < t1) t1 = t2;

	    bx += t1*ebx;
	    by += t1*eby;
	    my_state.bt += t1;
	    my_state.nx = cy - by;
	    my_state.ny = bx - cx;

	    continue;
	}

	if (tq < 0) {
	    //B stays fixed, C moves along ec
	    //possible stopping conditions: 
	    //   (1) C hits end of edge,
	    //   (2) BC becomes parallel to ea,
	    //   (3) tq becomes 0
	    t1 = 1 - my_state.ct;

	    t2 = -(by*eax - cy*eax - bx*eay + cx*eay)/(eay*ecx - eax*ecy);

	    t3  = ebx*ecx*(by-cy)*(by-cy);
	    t3 += eby*ecy*(bx-cx)*(bx-cx);
	    t3 += ebx*ecy*(-ay*bx + ax*by - bx*by + ay*cx - ax*cy + 2*bx*cy - cx*cy);
	    t3 += eby*ecx*( ay*bx - ax*by - bx*by - ay*cx + 2*by*cx + ax*cy - cx*cy);
	    t3 /= -ee*fc;

	    if (t3 < t1) t1 = t3;
	    if (t2 < t1) t1 = t2;

	    cx += t1*ecx;
	    cy += t1*ecy;
	    my_state.ct += t1;
	    my_state.nx = cy - by;
	    my_state.ny = bx - cx;

	    continue;
	}

	if (tq == 0) {
	    //B and C both move
	    //possible stopping conditions: 
	    //   (1) B hits end of edge,
	    //   (2) C hits end of edge,
	    //   (3) BC becomes parallel to ea,

	    //time in terms of b's motion until it hits vertex
	    t1 = 1 - my_state.bt;

	    //time in terms of c's motion until it hits vertex
	    t2 = 1 - my_state.ct;
	    //time in terms of b's motion until c hits vertex
	    t2 = (ay*ecx - 2*by*ecx + cy*ecx - ax*ecy + 2*bx*ecy - cx*ecy)*t2/(ay*ebx + by*ebx - 2*cy*ebx - ax*eby - bx*eby + 2*cx*eby - 2*ee*t2);

	    t3  = ebx*ecx*eay*(by-cy);
	    t3 += ebx*ecy*(-ay*eax + by*eax + ax*eay - 2*bx*eay + cx*eay);
	    t3 += eby*ecx*( ay*eax - 2*by*eax + cy*eax - ax*eay + bx*eay);
	    t3 += eby*ecy*eax*(bx-cx);

	    //time in terms of b's motion until (1 or 2 or 3)
	    if (t3 < t1) t1 = t3;
	    if (t2 < t1) t1 = t2;
	    //time in terms of c's motion until (1 or 2 or 3)
	    t2 = (ay*ebx + by*ebx - 2*cy*ebx - ax*eby - bx*eby + 2*cx*eby)*t1/(ay*ecx - 2*by*ecx + cy*ecx - ax*ecy + 2*bx*ecy - cx*ecy + 2*ee*t1);

	    bx += t1*ebx;
	    by += t1*eby;
	    my_state.bt += t1;
	    cx += t2*ecx;
	    cy += t2*ecy;
	    my_state.ct += t2;
	    my_state.nx = cy - by;
	    my_state.ny = bx - cx;
	}
    }

    my_state.ai = ret[0];
    my_state.bi = ret[1];
    my_state.ci = ret[2];
    return my_state;


}
