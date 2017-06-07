#include <string>
#include <fstream>
#include "max-triangle.hpp"

bool test_anchor(){
    std::vector <std::vector <mpq_class> > polygons = {{3040, 4460, 2506, 4423, 759, 2927, 1000, 1000, 1213, 691, 3383, 413, 5000, 1000, 4752, 4262, 4745, 4322},{707, 707, 258, 965, -258, 965, -707, 707, -965, 258, -965, -258, -707, -707, -258, -965, 258, -965, 707, -707, 965, -258, 965, 258}};
    std::vector <mpq_class> normals = {1000, 0, 965, 258, 866, 500, 707, 707, 500, 866, 258, 965, 0, 1000, -258, 965, -500, 866, -707, 707, -866, 500, -965, 258, -1000, 0, -965, -258, -866, -500, -707, -707, -500, -866, -258, -965, 0, -1000, 258, -965, 500, -866, 707, -707, 866, -500, 965, -258};
    std::vector <std::string> answers = { "4745/1", "489105/539", "4745/1", "4322/1", "5000/1", "1000/1", "6596369440/1609721", "7041532300/1609721", "5000/1", "1000/1", "90290000/30059", "133987660/30059", "22731048614502/4573535527", "12741295770649/9147071054", "14442431582001/6561369482", "13653989196599/3280684741", "10042817138245824/2049207949003", "4722240784813694/2049207949003", "2186256639044071/1097120460666", "2185171897579159/548560230333", "22709194554162112/4679521977703", "13734404319070422/4679521977703", "10004525165592089/5749137493126",
       	"10829097693494121/2874568746563", "23877745938624/4949858291", "10063453874/3034861", "5508161685573/4540152056", "10063453874/3034861", "7663272108/1605907", "6423445523/1605907", "759/1", "2927/1", "4745/1", "4322/1", "773616584/894641", "1862334152/894641", "25762434667740/5993946511", "26122658011580/5993946511", "875478988740/881367067", "928447261220/881367067", "166253000/55493", "247329320/55493", "1000/1", "1000/1", "2506/1", "4423/1", "1622646496/1082887", "708675346/1082887", "6618173049/3731986", "12373322924617/3259889771", "6618173049/3731986",
       	"358187633343/578457830", "5514173436902064/3866007215239", "13525014020527937/3866007215239", "4171955270136522/1886823079859", "1062531887792532/1886823079859", "378482230637456/301344807549", "1010280806060443/301344807549", "473911704773562/162361916923", "383548573951044/811809584615", "2653607/3243", "9656821/3243", "3383/1", "413/1", "759/1", "2927/1", "2926890997/658421", "525841767/658421", "1548665000/1797377", "3786041000/1797377", "5000/1", "1000/1", "1000/1", "1000/1", "5000/1", "1000/1", "1213/1", "691/1", "7945356236/1605907", "2713129291/1605907",
       	"1213/1", "691/1", "3585308828/737223", "2063149343/737223", "14950948687/7687944", "142342285753/238326264", "5894429636/1224795", "131565767171/37968645", "85764515/33637", "87414757/168185", "4745/1", "4322/1", "584555537/156551", "253803121/469653", "4745/1", "4322/1", "250879/516", "-748819/898", "250879/516", "748819/898", "352536102/499849", "-353885765/499849", "258/1", "965/1", "41969661991/50279918", "-108561463/223964", "-448096/194417", "965/1", "965/1", "-258/1", "-258/1", "965/1", "965/1", "-448096/194417", "-108561463/223964", "41969661991/50279918",
       	"965/1", "258/1", "-353885765/499849", "352536102/499849", "748819/898", "250879/516", "-748819/898", "250879/516", "353885765/499849", "352536102/499849", "-965/1", "258/1", "108561463/223964", "41969661991/50279918", "-965/1", "-448096/194417", "258/1", "965/1", "-965/1", "-258/1", "448096/194417", "965/1", "-41969661991/50279918", "-108561463/223964", "-258/1", "965/1", "-352536102/499849", "-353885765/499849", "-250879/516", "748819/898", "-250879/516", "-748819/898", "-352536102/499849", "353885765/499849", "-258/1", "-965/1", "-41969661991/50279918",
       	"108561463/223964", "448096/194417", "-965/1", "-965/1", "258/1", "258/1", "-965/1", "-965/1", "448096/194417", "108561463/223964", "-41969661991/50279918", "-965/1", "-258/1", "353885765/499849", "-352536102/499849", "-748819/898", "-250879/516", "748819/898", "-250879/516", "-353885765/499849", "-352536102/499849", "965/1", "-258/1", "-108561463/223964", "-41969661991/50279918", "965/1", "448096/194417", "-258/1", "-965/1", "965/1", "258/1", "-448096/194417", "-965/1", "41969661991/50279918", "108561463/223964", "258/1", "-965/1", "352536102/499849", "353885765/499849"};

    unsigned int i,j,m;
    struct state my_state;

    m=normals.size()/2;

    for (i = 0; i < polygons.size(); i++) {
	for (j = 0; j < m; j++) {
	    my_state = anchored_triangle(polygons[i], normals[2*j+0], normals[2*j+1]);
	    my_state.ai += polygons[i].size()/2;
	    my_state.ai %= polygons[i].size()/2;
	    my_state.bi += polygons[i].size()/2;
	    my_state.bi %= polygons[i].size()/2;
	    my_state.ci += polygons[i].size()/2;
	    my_state.ci %= polygons[i].size()/2;
	    if (my_state.status != status_ok) return false;
	    if (mpq_class(answers[i*m*4 + j*4 + 0]) != polygons[i][2*my_state.bi+0]*(1-my_state.bt)+polygons[i][(2*my_state.bi+2)%polygons[i].size()]*my_state.bt) return false;
	    if (mpq_class(answers[i*m*4 + j*4 + 1]) != polygons[i][2*my_state.bi+1]*(1-my_state.bt)+polygons[i][(2*my_state.bi+3)%polygons[i].size()]*my_state.bt) return false;
	    if (mpq_class(answers[i*m*4 + j*4 + 2]) != polygons[i][2*my_state.ci+0]*(1-my_state.ct)+polygons[i][(2*my_state.ci+2)%polygons[i].size()]*my_state.ct) return false;
	    if (mpq_class(answers[i*m*4 + j*4 + 3]) != polygons[i][2*my_state.ci+1]*(1-my_state.ct)+polygons[i][(2*my_state.ci+3)%polygons[i].size()]*my_state.ct) return false;
	    //std::cout << "A: " << polygons[i][2*my_state.ai+0] << ", " << polygons[i][2*my_state.ai+1] << std::endl;
	    //std::cout << "B: " << polygons[i][2*my_state.bi+0]*(1-my_state.bt)+polygons[i][(2*my_state.bi+2)%polygons[i].size()]*my_state.bt  << ", " << polygons[i][2*my_state.bi+1]*(1-my_state.bt)+polygons[i][(2*my_state.bi+3)%polygons[i].size()]*my_state.bt << std::endl;
	    //std::cout << "C: " << polygons[i][2*my_state.ci+0]*(1-my_state.ct)+polygons[i][(2*my_state.ci+2)%polygons[i].size()]*my_state.ct  << ", " << polygons[i][2*my_state.ci+1]*(1-my_state.ct)+polygons[i][(2*my_state.ci+3)%polygons[i].size()]*my_state.ct << std::endl;
 	}
    }

    return true;
}

bool test_max() {
    std::vector <std::vector <mpq_class> > polygons;

    unsigned int i;
    struct state s1,s2;
    mpq_class ax,ay,bx,by,cx,cy;
    mpq_class z1,z2;

    std::string line, entry;
    std::ifstream infile("test-polygons");

    while (std::getline(infile, line)) {
	std::stringstream line_stream(line);
	polygons.push_back(std::vector <mpq_class>());
	while (line_stream >> entry) {
	    if (!entry.empty()) polygons.back().push_back(mpq_class(entry));
	}
    }

    for (i = 0; i < polygons.size(); i++) {
    //for (i = 28; i < 29; i++) {
	s1 = maximum_triangle(polygons[i]);
	if (s1.status != status_ok) {
	    std::cout << i << " ";
	    print_status(s1.status);
	    return false;
	}
	s2 = naive_maximum_triangle(polygons[i]);
	ax = polygons[i][2*s1.ai+0];
	ay = polygons[i][2*s1.ai+1];
	bx = polygons[i][2*s1.bi+0];
	by = polygons[i][2*s1.bi+1];
	cx = polygons[i][2*s1.ci+0];
	cy = polygons[i][2*s1.ci+1];
	z1 = (bx-ax)*(cy-ay) - (cx-ax)*(by-ay);
    	ax = polygons[i][2*s2.ai+0];
	ay = polygons[i][2*s2.ai+1];
	bx = polygons[i][2*s2.bi+0];
	by = polygons[i][2*s2.bi+1];
	cx = polygons[i][2*s2.ci+0];
	cy = polygons[i][2*s2.ci+1];
	z2 = (bx-ax)*(cy-ay) - (cx-ax)*(by-ay);
	std::cout << z1 << " = " << z2 << std::endl;
	if (z1 != z2) return false;
    }
    return true;
}

int main(int argc, char *argv[]){
    if (test_anchor()) std::cout << "anchor triangle finding unit test PASSED" << std::endl;
    else std::cout << "anchor triangle finding unit test FAILED" << std::endl;
    if (test_max()) std::cout << "maximum triangle finding unit test PASSED" << std::endl;
    else std::cout << "maximum triangle finding unit test FAILED" << std::endl;
    
    /*std::vector <std::vector <mpq_class> > polygons = {{3040, 4460, 2506, 4423, 759, 2927, 1000, 1000, 1213, 691, 3383, 413, 5000, 1000, 4752, 4262, 4745, 4322},{707, 707, 258, 965, -258, 965, -707, 707, -965, 258, -965, -258, -707, -707, -258, -965, 258, -965, 707, -707, 965, -258, 965, 258}};

    struct state my_state = maximum_triangle(polygons[0]);
    if (my_state.status != status_ok) std::cout << "status not OK" << std::endl;
    else std::cout << my_state.ai << "\t" << my_state.bi << "\t" << my_state.ci << std::endl;*/
    return 0;
}

