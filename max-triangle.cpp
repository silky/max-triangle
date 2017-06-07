#include <iomanip>
#include "max-triangle.hpp"

#define getmod(arg1,arg2) arg1[(((arg2)+(arg1.size()))%(arg1.size()))]

void print_status(state_flag status) {
    switch(status) {
	case status_ok: std::cout << "status_ok" << std::endl; break;
	case status_not_convex: std::cout << "status_not_convex" << std::endl; break;
	case status_no_interior: std::cout << "status_no_interior" << std::endl; break;
	case status_runtime_error: std::cout << "status_runtime_error" << std::endl; break;
	case status_max_not_at_vertex: std::cout << "status_max_not_at_vertex" << std::endl; break;
	case status_maxiter_exceeded: std::cout << "status_maxiter_exceeded" << std::endl; break;
    }
}

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

    struct state my_state;
    int ret[3] = {0,0,0};
    unsigned int sz = polygon.size()/2;
    bool max_at_vert=false;
    unsigned int iter = 0;
    unsigned int maxiter = sz*10;
    int ai_start;

    mpq_class ax,ay,bx,by,cx,cy;
    mpq_class pax,pay,pbx,pby,pcx,pcy;
    mpq_class eax,eay,ebx,eby,ecx,ecy;
    mpq_class t1,t2,t3;
    mpq_class ee,fb,fc,tq;
    mpq_class area;

    my_state = anchored_triangle(polygon,polygon[1]-polygon[3],polygon[2]-polygon[0]);
    my_state.nx = polygon[1]-polygon[3];
    my_state.ny = polygon[2]-polygon[0];
    if (my_state.status != status_ok) {return my_state;}
    ai_start = my_state.ai = 0;

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

    bx += my_state.bt*ebx;
    by += my_state.bt*eby;
    cx += my_state.ct*ecx;
    cy += my_state.ct*ecy;

    //mpq_class amax = (bx-ax)*(cy-ay) - (cx-ax)*(by-ay);
    //ret[0] = my_state.ai; ret[1] = my_state.bi; ret[2] = my_state.ci;
    mpq_class amax = -1;


    while (my_state.ai <= ai_start + (int) sz) {

	if (iter++ > maxiter) {my_state.status=status_maxiter_exceeded; return my_state;}

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
	
#ifdef DEBUG
	std::cout << my_state.ai << "\t" << my_state.bi << "\t" << my_state.ci << "\t";
        std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << my_state.bt.get_d() << "\t";
	std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << my_state.ct.get_d() << "\t";
	std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << area.get_d() << std::endl;
#endif


	ee = ebx*ecy - ecx*eby;
	fb = ebx*(ay + by - 2*cy) - eby*(ax + bx - 2*cx);
	fc = ecx*(ay + cy - 2*by) - ecy*(ax + cx - 2*bx);
	tq  = my_state.nx*(by*ebx*ecx - cy*ebx*ecx - ax*eby*ecx + cx*eby*ecx + ax*ebx*ecy - bx*ebx*ecy);
	tq -= my_state.ny*(ay*eby*ecx - by*eby*ecx - ay*ebx*ecy + cy*ebx*ecy + bx*eby*ecy - cx*eby*ecy);
#ifdef DEBUG
	std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << tq.get_d() << "\t";
	std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << ee.get_d() << "\t";
	std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << fb.get_d() << "\t";
	std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << fc.get_d() << std::endl;
#endif
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

	    if (eay*ebx - eax*eby == 0) t2 = -1;
	    else t2 = (by*eax - cy*eax - bx*eay + cx*eay)/(eay*ebx - eax*eby);

	    if (ee == 0) t3 = -1;
	    else {
		t3  = ebx*ecx*(by-cy)*(by-cy);
		t3 += eby*ecy*(bx-cx)*(bx-cx);
		t3 += ebx*ecy*(-ay*bx + ax*by - bx*by + ay*cx - ax*cy + 2*bx*cy - cx*cy);
		t3 += eby*ecx*( ay*bx - ax*by - bx*by - ay*cx + 2*by*cx + ax*cy - cx*cy);
		t3 /= ee*fb;
	    }

#ifdef DEBUG
	    std::cout << "tq: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << tq.get_d() << "\t";
	    std::cout << "t1: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t1.get_d() << "\t";
	    std::cout << "t2: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t2.get_d() << "\t";
	    std::cout << "t3: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t3.get_d() << std::endl;
#endif

	    if (t3 < t1 && t3 > 0) t1 = t3;
	    if (t2 < t1 && t2 > 0) t1 = t2;

	    if (t1 < 0) {my_state.status=status_runtime_error; return my_state;}

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

	    //std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << ecy.get_d() << "\t";
	    //std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << eax.get_d() << "\t";
	    //std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << eay.get_d() << std::endl;
	    //std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << cx.get_d() << "\t";
	    //std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << cy.get_d() << "\t";
	    //std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << bx.get_d() << "\t";
	    //std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << by.get_d() << std::endl;
	    if (eay*ecx - eax*ecy == 0) t2 = -1;
	    else t2 = -(by*eax - cy*eax - bx*eay + cx*eay)/(eay*ecx - eax*ecy);

	    if (ee == 0) t3 = -1;
	    else {
		t3  = ebx*ecx*(by-cy)*(by-cy);
		t3 += eby*ecy*(bx-cx)*(bx-cx);
		t3 += ebx*ecy*(-ay*bx + ax*by - bx*by + ay*cx - ax*cy + 2*bx*cy - cx*cy);
		t3 += eby*ecx*( ay*bx - ax*by - bx*by - ay*cx + 2*by*cx + ax*cy - cx*cy);
		t3 /= -ee*fc;
	    }

#ifdef DEBUG
	    std::cout << "tq: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << tq.get_d() << "\t";
	    std::cout << "t1: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t1.get_d() << "\t";
	    std::cout << "t2: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t2.get_d() << "\t";
	    std::cout << "t3: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t3.get_d() << std::endl;
#endif


	    if (t3 < t1 && t3 > 0) t1 = t3;
	    if (t2 < t1 && t2 > 0) t1 = t2;
	    
	    if (t1 < 0) {my_state.status=status_runtime_error; return my_state;}

	    cx += t1*ecx;
	    cy += t1*ecy;
	    my_state.ct += t1;
	    my_state.nx = cy - by;
	    my_state.ny = bx - cx;

	    continue;
	}

	if (tq == 0) {
	    mpq_class t1b,t1c,t2b,t2c,t3b,t3c;
	    //B and C both move
	    //possible stopping conditions: 
	    //   (1) B hits end of edge,
	    //   (2) C hits end of edge,
	    //   (3) BC becomes parallel to ea,
#ifdef DEBUG
	    mpq_class qu = (ay*ebx + by*ebx - 2*cy*ebx - ax*eby - bx*eby + 2*cx*eby);
	    std::cout << "alpha1: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << qu.get_d() << "\t";
	    qu = (ay*ecx - 2*by*ecx + cy*ecx - ax*ecy + 2*bx*ecy - cx*ecy);
	    std::cout << "alpha2: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << qu.get_d() << "\t";
	    qu = (eay*ebx - eax*eby);
	    std::cout << "eaxab: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << qu.get_d() << "\t";
	    std::cout << "ee: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << ee.get_d() << std::endl;
#endif

	    //time in terms of b's motion until it hits vertex
	    t1b = 1 - my_state.bt;
	    if (ay*ecx - 2*by*ecx + cy*ecx - ax*ecy + 2*bx*ecy - cx*ecy + 2*ee*t1b == 0) {t1b = t1c = -1;}
	    else t1c = (ay*ebx + by*ebx - 2*cy*ebx - ax*eby - bx*eby + 2*cx*eby)*t1b/(ay*ecx - 2*by*ecx + cy*ecx - ax*ecy + 2*bx*ecy - cx*ecy + 2*ee*t1b);

	    //time in terms of c's motion until it hits vertex
	    t2c = 1 - my_state.ct;
	    //time in terms of b's motion until c hits vertex
	    if (ay*ebx + by*ebx - 2*cy*ebx - ax*eby - bx*eby + 2*cx*eby - 2*ee*t2c == 0) {t2b = t2c = -1;}
	    else t2b = (ay*ecx - 2*by*ecx + cy*ecx - ax*ecy + 2*bx*ecy - cx*ecy)*t2c/(ay*ebx + by*ebx - 2*cy*ebx - ax*eby - bx*eby + 2*cx*eby - 2*ee*t2c);

	    if (eay*ebx - eax*eby == 0) {t3b = t3c = -1;}
	    else if (eay*ecx - eax*ecy == 0) {t3b = t3c = -1;}
	    else {
		t3b  = ebx*ecx*eay*(by-cy);
		t3b += ebx*ecy*(-ay*eax + by*eax + ax*eay - 2*bx*eay + cx*eay);
		t3b += eby*ecx*( ay*eax - 2*by*eax + cy*eax - ax*eay + bx*eay);
		t3b += eby*ecy*eax*(bx-cx);
		t3b /= 2*ee*(eay*ebx - eax*eby);

		t3c  = ebx*ecx*eay*(by-cy);
		t3c += ebx*ecy*(-ay*eax - by*eax + 2*cy*eax + ax*eay - cx*eay);
		t3c += eby*ecx*( ay*eax - cy*eax - ax*eay - bx*eay + 2*cx*eay);
		t3c += eby*ecy*eax*(bx-cx);
		t3c /= 2*ee*(eay*ecx - eax*ecy);
	    }

#ifdef DEBUG
	    std::cout << "tq: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << tq.get_d() << "\t";
	    std::cout << "t1b: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t1b.get_d() << "\t";
	    std::cout << "t2b: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t2b.get_d() << "\t";
	    std::cout << "t3b: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t3b.get_d() << std::endl;
	    std::cout << "tq: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << tq.get_d() << "\t";
	    std::cout << "t1c: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t1c.get_d() << "\t";
	    std::cout << "t2c: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t2c.get_d() << "\t";
	    std::cout << "t3c: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t3c.get_d() << std::endl;
#endif

	    //time in terms of b's motion until (1 or 2 or 3)
	    /*if (t2 < t1) t1 = t2;
	    if (t3 < t1 && t3 > 0) t1 = t3;
	    if (t1 < 0) {my_state.status=status_runtime_error; return my_state;}*/
	    
	    //if (t2b <= t1b && t2c <= t1c) {t1b = t2b; t1c = t2c;}
	    if ( (t1c < 0) || (t2b <= t1b && t2c <= t1c) ) { t1c = t2c; t1b = t2b;}
	    if (t3b <= t1b && t3c <= t1c && t3b >= 0 && t3c >= 0) {t1b = t3b; t1c = t3c;}
	    if (t1b < 0 || t1c < 0) {my_state.status=status_runtime_error; return my_state;}
	    //time in terms of c's motion until (1 or 2 or 3)
	    //t2 = (ay*ebx + by*ebx - 2*cy*ebx - ax*eby - bx*eby + 2*cx*eby)*t1/(ay*ecx - 2*by*ecx + cy*ecx - ax*ecy + 2*bx*ecy - cx*ecy + 2*ee*t1);

	    bx += t1b*ebx;
	    by += t1b*eby;
	    my_state.bt += t1b;
	    cx += t1c*ecx;
	    cy += t1c*ecy;
	    my_state.ct += t1c;
	    my_state.nx = cy - by;
	    my_state.ny = bx - cx;
	}
    }

    if (!max_at_vert) my_state.status = status_max_not_at_vertex;
    my_state.ai = ret[0];
    my_state.bi = ret[1];
    my_state.ci = ret[2];
    return my_state;


}

struct state naive_maximum_triangle(std::vector<mpq_class> polygon){
    mpq_class amax, area;
    mpq_class ax,ay,bx,by,cx,cy;
    unsigned int sz = polygon.size()/2;
    int ret[3] = {0,0,0};
    struct state my_state;

    amax = -1;
    for (my_state.ai = 0; my_state.ai+2 < (int) sz; my_state.ai++) { 
	ax = getmod(polygon, 2*my_state.ai+0);
	ay = getmod(polygon, 2*my_state.ai+1);
	for (my_state.bi = my_state.ai+1; my_state.bi+1 < (int) sz; my_state.bi++) { 
	    bx = getmod(polygon, 2*my_state.bi+0);
	    by = getmod(polygon, 2*my_state.bi+1);
	    for (my_state.ci = my_state.bi+1; my_state.ci < (int) sz; my_state.ci++) { 
		cx = getmod(polygon, 2*my_state.ci+0);
		cy = getmod(polygon, 2*my_state.ci+1);
		area = (bx-ax)*(cy-ay) - (cx-ax)*(by-ay);
		if (area > amax) {
		    amax = area;
		    ret[0] = my_state.ai;
		    ret[1] = my_state.bi;
		    ret[2] = my_state.ci;
		}
	    }
	}
    }
    my_state.ai = ret[0];
    my_state.bi = ret[1];
    my_state.ci = ret[2];
    return my_state;
}
