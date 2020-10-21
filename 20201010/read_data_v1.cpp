#include <iostream>
#include <string>
#include <bitset>
#include <chrono>
#include <fstream>
#include <cmath>

using namespace std;

class monitor_constant {
/*
* This class have the constant array of two 16PU. These array will be changed to Eigen::Matrix.
* Basically, 16PU outputs twiss parameter using formula, which is,
*               Voltage = Matrix * Momentum
* , so to obtain result, it need those monitor parameter.
*   Matrix(inversed) 16x16 : Calculated with Wire Calibratoin. Using this, we can get beam parameter such as positon, or emittance.
*   Gain 16 : For fixing voltage gain error
*   Align : For fixing position error
* ISSUE :
* 1) Because of lack of skill of C++, when I tried to compile this code, it return several warning of in-class initialization.
*    But compiled code run correctly(maybe), I didn't touch anything about this.
*/
    public:
        const double matrix13[16][16] = 
        {{9.999413998070250109e-01,9.999460435510404421e-01,1.000036939761390409e+00,9.998124202563855034e-01,1.000271488239349749e+00,9.998693920504219124e-01,1.000111480706963318e+00,9.998456629174904409e-01,9.999773275920876836e-01,1.000162081686137849e+00,9.998997463053272972e-01,1.000186489373741328e+00,9.997085318373660767e-01,1.000132905206968958e+00,9.999092115071314124e-01,1.000191401767124733e+00},
        {8.402690869148989350e+01,7.855569013495049546e+01,5.868037260691073698e+01,3.282522159378924442e+01,-7.166513112806041086e-01,-3.193109828147644436e+01,-6.147302662371219384e+01,-7.772471435037130050e+01,-8.545037779220223229e+01,-7.759246361210819032e+01,-6.062931303009958839e+01,-3.158140091037386199e+01,-1.446554977717688650e-01,3.261407689753472994e+01,5.948528506596413479e+01,7.875739354011824389e+01},
        {-3.184580464191904547e-01,3.228656153914549520e+01,5.995639198359812383e+01,7.820350744773398333e+01,8.476220176545791674e+01,7.757651271415326732e+01,6.031722793075426381e+01,3.220030029154117557e+01,4.337302640585219415e-01,-3.240834556455450866e+01,-5.952089441992973917e+01,-7.802962192320494239e+01,-8.487089293163909076e+01,-7.812708433428798571e+01,-5.990049734393068093e+01,-3.223788653676192695e+01},
        {7.137783296339693152e+03,5.112365477098272095e+03,-7.851964983994369618e+01,-5.081306374642089395e+03,-7.195729819149987634e+03,-5.035351301641623650e+03,5.802079164518411147e+01,5.085147097495299931e+03,7.223963968785255929e+03,5.052380434869687633e+03,3.093289319604604870e+01,-5.086032890206775846e+03,-7.219966991685299035e+03,-5.069979210138622875e+03,-4.370857235080535474e+00,5.135855289893223926e+03},
        {-4.128864486652737753e+01,5.087607693875425866e+03,7.136140344759595791e+03,5.099855744849646726e+03,-4.102723224737526664e+01,-5.040457598271681491e+03,-7.316189874531491114e+03,-5.057047774689844118e+03,-2.612702648095267577e+01,5.048683556148226671e+03,7.205862535473165735e+03,5.006432193559158804e+03,-2.916023524624746344e+01,-5.076237364001140122e+03,-7.202681364551679508e+03,-5.072362732333781423e+03},
        {6.049011323557938449e+05,2.396483519005203561e+05,-4.382757016334204818e+05,-5.618301118111484684e+05,-1.571629856130904727e+03,5.632854633836707799e+05,4.344627039328728570e+05,-2.369006016064683208e+05,-6.115075652575897984e+05,-2.343046354074297997e+05,4.341596774890003726e+05,5.562759611219625222e+05,-2.106216762353617469e+03,-5.657101597251122585e+05,-4.308280296919397661e+05,2.393924747899266367e+05},
        {-2.662443326703010825e+03,5.646673200289527886e+05,4.279343994867376168e+05,-2.366162533193831914e+05,-6.090826484717677813e+05,-2.330150638983670797e+05,4.439187612618173589e+05,5.646857144498244161e+05,-2.418693814792285139e+03,-5.565914658578806557e+05,-4.390198956733581144e+05,2.405489163061128347e+05,6.118203734236392193e+05,2.341006563270064944e+05,-4.343683409772803425e+05,-5.654667673793170834e+05},
        {5.087124571121729910e+07,1.119067933127190452e+06,-5.272955327068960667e+07,8.680483496539668413e+05,5.111064435668140650e+07,5.821473845092523843e+05,-5.345720835225655884e+07,1.006413324495706591e+06,5.141201192359492183e+07,7.237471585985765560e+05,-5.296820200364012271e+07,1.323752426483716350e+06,5.158727877359198779e+07,2.446428717375808046e+05,-5.230182401809504628e+07,8.855952866770752007e+05},
        {-3.838756986744244932e+05,5.241258514119739085e+07,-6.772074381322815316e+05,-5.168932041957750171e+07,-5.462655290799565846e+05,5.222136688028974086e+07,-8.780164960629242705e+05,-5.219173148805026710e+07,5.996072744596879929e+05,5.105952124851896614e+07,4.918003338072042679e+05,-5.171546453306925297e+07,3.770671803330594557e+05,5.209912062451072782e+07,-4.198797444921328570e+05,-5.208505744134169817e+07},
        {4.304301365307842255e+09,-1.599427431200844288e+09,-3.150392094739248276e+09,4.117372586963991642e+09,1.440313067957474105e+07,-4.085007020188084126e+09,3.254280087345081806e+09,1.586597720981162310e+09,-4.327028667759906769e+09,1.603631306629714966e+09,3.173473113201637268e+09,-4.082679482509153366e+09,3.369625028280577809e+07,4.095570413339972973e+09,-3.184360063912318707e+09,-1.595812367223905325e+09},
        {-2.102333120259914547e+07,4.125340881569543839e+09,-3.180174857363269329e+09,-1.615202780518994808e+09,4.320756676590899467e+09,-1.623865802434098244e+09,-3.183257980723783493e+09,4.136961603416942120e+09,-3.777085914116869122e+07,-4.054970198830155849e+09,3.158694890463657379e+09,1.568760697542835236e+09,-4.363833028166501999e+09,1.656159573613260984e+09,3.124283466705462456e+09,-4.116477574751970768e+09},
        {3.604535387839906006e+11,-2.554800833334613037e+11,-1.207473219268465191e+08,2.585142133768315430e+11,-3.620215949421724854e+11,2.553867624233882751e+11,-2.375439299573695183e+09,-2.564574463999008484e+11,3.632380267932874756e+11,-2.560974503151954346e+11,1.396151675565687180e+09,2.512340434692117004e+11,-3.650361730757618408e+11,2.594215007235646667e+11,-3.893049071341799259e+09,-2.551705562390465393e+11},
        {2.143926967192851603e+08,2.665171992605645752e+11,-3.731634033223012085e+11,2.675063575278382263e+11,-6.247197764797392488e+07,-2.626556799043166199e+11,3.794801826879353638e+11,-2.692885313205568848e+11,2.790921686277306080e+09,2.623359363631508484e+11,-3.718211762495717773e+11,2.636847362527774353e+11,-1.689168636435148001e+09,-2.652555795022042542e+11,3.749221298006658325e+11,-2.695171217885629578e+11},
        {2.937521841317168750e+13,-2.745575140157400000e+13,2.125136478914786328e+13,-1.207462281257746680e+13,5.157480427142016602e+10,1.174723629863906055e+13,-2.147299785263950781e+13,2.753898628755692188e+13,-2.945194560154652734e+13,2.721046700469135938e+13,-2.110487736917602344e+13,1.182805539971532422e+13,-7.890405831542990112e+10,-1.188735452162009375e+13,2.120539213286899609e+13,-2.754674298073397266e+13},
        {5.591433517171932220e+10,1.195678964676207422e+13,-2.114999901277765234e+13,2.748473708065037109e+13,-2.919485683621238672e+13,2.703571649562666016e+13,-2.157422929283521875e+13,1.216316281476184961e+13,-1.794486419635328064e+11,-1.169456655227083203e+13,2.090942587083592188e+13,-2.675360293263904297e+13,2.942612648068332422e+13,-2.744143813504088672e+13,2.142687804228785156e+13,-1.221366315172071094e+13},
        {1.411664770215943750e+15,-1.425581589371573500e+15,1.401387927769147000e+15,-1.424603075057991750e+15,1.396734144148960750e+15,-1.395123838367636750e+15,1.417306274271324000e+15,-1.428928753756292000e+15,1.410059836059073250e+15,-1.402477989716191000e+15,1.376925171344775250e+15,-1.376245229453973500e+15,1.400404994852157000e+15,-1.414332897004222250e+15,1.406113595484454750e+15,-1.432234529429722500e+15}};
        const double matrix15[16][16] = 
        {{1.000003842111640173e+00,9.999339155758250053e-01,9.999613485487530706e-01,9.999245728760640661e-01,1.000169133027317470e+00,9.999887853817465144e-01,9.999875367206880750e-01,9.999073921609406002e-01,9.999856232906959574e-01,1.000073356874412234e+00,1.000017455261239752e+00,1.000071559698644696e+00,9.998422062121793408e-01,1.000044355136718854e+00,9.999446488560410629e-01,1.000140151052085535e+00},
        {8.380796005233882795e+01,7.771145478628254466e+01,5.962960636531391856e+01,3.198826049418720885e+01,4.229671786544587087e-01,-3.218673900362998808e+01,-5.958301014935136664e+01,-7.845578013674658280e+01,-8.456861448758583322e+01,-7.793976855751184019e+01,-6.089376232309970760e+01,-3.217744786098418075e+01,-3.785840078649872070e-02,3.238496513906805774e+01,6.003392960722874960e+01,7.779868772824023893e+01},
        {9.622955467608648172e-01,3.199108449008893373e+01,6.197620014722254211e+01,7.764073970663888247e+01,8.666517794977795575e+01,7.749201832231049991e+01,6.112833626807520915e+01,3.181451876045507277e+01,7.090234546601680288e-01,-3.312628849943363463e+01,-5.958684108496061782e+01,-7.806307359498298126e+01,-8.412913589212907084e+01,-7.821944695427288252e+01,-5.922829736390616517e+01,-3.306211160653684544e+01},
        {7.018994189263068620e+03,5.080357160639074209e+03,-1.837478689975672523e+02,-5.080218438777588744e+03,-7.366222460595396115e+03,-5.062213150046626652e+03,-6.330767817191470925e+01,5.136716320112727772e+03,7.194531204534600874e+03,5.060009952214029909e+03,4.143149338333218168e+01,-5.025487138487653283e+03,-7.214379876305093603e+03,-5.014752753263587692e+03,-2.908876472634371169e+01,5.075261854226389005e+03},
        {4.719073310244203157e+01,5.068121272655576831e+03,7.272526950383818985e+03,5.026391884011711227e+03,3.344636575995555461e+01,-5.056848282609405942e+03,-7.228762161704264145e+03,-5.107356320504989526e+03,1.881793136202152184e+01,5.064331720763235353e+03,7.364170437630719789e+03,5.003562853610916136e+03,2.785927761432990835e+01,-5.070438958997464397e+03,-7.198569966895035577e+03,-5.081946293932432127e+03},
        {5.931658195265893592e+05,2.338815616208376305e+05,-4.476826783846946782e+05,-5.602819832985377871e+05,-1.639226467128859667e+03,5.638122697122534737e+05,4.375989044537804439e+05,-2.385500737585728930e+05,-6.142313662606854923e+05,-2.341003643493511190e+05,4.422583155467043398e+05,5.564866324803766329e+05,1.553429189010804066e+03,-5.593735328294278588e+05,-4.361454774954200839e+05,2.352126008516644652e+05},
        {1.733590102775183823e+03,5.636934467244467232e+05,4.260196134581808583e+05,-2.405033210676578165e+05,-6.225448630826633889e+05,-2.370911855171922653e+05,4.360538632523380802e+05,5.695775586681684945e+05,-3.431410642116015424e+03,-5.617377246291603660e+05,-4.498123185435839114e+05,2.349160907993594592e+05,6.111627218046365306e+05,2.317092636066973209e+05,-4.345437074614539160e+05,-5.605178167293900624e+05},
        {4.987092356801528484e+07,4.447518450176328188e+05,-5.295405933303615451e+07,1.458359062735761749e+06,5.220777192180593312e+07,9.405731605342011899e+05,-5.344339141771052033e+07,9.279005807395989541e+05,5.188593210407877713e+07,6.072565646811689949e+05,-5.424851930740482360e+07,8.082962582123586908e+05,5.162304938693358004e+07,3.637025648001093650e+05,-5.258282670158503950e+07,6.282779393843616126e+05},
        {1.977806225841316627e+05,5.169505295087503642e+07,-1.130016897219566396e+06,-5.212425412875749916e+07,1.357198185928351886e+05,5.213947339811725169e+07,-8.250188782278122380e+04,-5.270661495766038448e+07,2.963350411871391698e+05,5.208375103712174296e+07,2.079337748516576830e+05,-5.132615877065725625e+07,-2.330371215148026240e+05,5.161127541797751188e+07,-8.508464904555617250e+04,-5.154316436626729369e+07},
        {4.202112711203301430e+09,-1.634195666702368498e+09,-3.106909989698691845e+09,4.149427056621729374e+09,-2.202770379153316468e+07,-4.123984889188008308e+09,3.228289590372787476e+09,1.621709320575987101e+09,-4.400534968101716042e+09,1.639161894067622900e+09,3.265042470884566307e+09,-4.060240162135756016e+09,-8.020696722210187465e+06,4.041783743209918499e+09,-3.155110322165650368e+09,-1.596973679787877798e+09},
        {3.388537066449384391e+07,4.030980185253132343e+09,-3.205676374572342873e+09,-1.580978307620392799e+09,4.430568934738506317e+09,-1.631733733018624544e+09,-3.193446903317361832e+09,4.156792061131393433e+09,-4.499298766080981120e+06,-4.133607259863242626e+09,3.253068793993239403e+09,1.589668265389823675e+09,-4.350170441346466064e+09,1.605656955020836115e+09,3.176162853574741364e+09,-4.074707600644390583e+09},
        {3.491852905899704590e+11,-2.530942426647469788e+11,3.492940201083008766e+09,2.569342206703104248e+11,-3.726711443960348511e+11,2.604355857435857544e+11,-3.275830177691141605e+09,-2.585089495585712585e+11,3.693232092650638428e+11,-2.597955839410419617e+11,-6.870804176606903076e+08,2.535960551532515564e+11,-3.641042681879528198e+11,2.526051257965307312e+11,6.335898619488921165e+08,-2.521209547659824524e+11},
        {3.420532074989078522e+09,2.576447819911822815e+11,-3.709081933430524902e+11,2.707143839164924622e+11,-2.477324968804492950e+09,-2.664755235872942810e+11,3.794037234303944092e+11,-2.698843546482242737e+11,-6.309554232246243395e+06,2.685989224207263184e+11,-3.850580728768602905e+11,2.633616879409482117e+11,5.472732178109289408e+08,-2.612775108673136902e+11,3.735300257592744141e+11,-2.641227410725647278e+11},
        {2.793837643239760938e+13,-2.656053689846369531e+13,2.106818956430916797e+13,-1.220541616123223242e+13,1.410053179974940491e+11,1.192487904183614062e+13,-2.146026365579735938e+13,2.771347116307352734e+13,-3.001522284013966016e+13,2.777374808534591797e+13,-2.191481854220090625e+13,1.191792536429475586e+13,-8.449977767506344604e+10,-1.155660444583321289e+13,2.093750252803382422e+13,-2.663338881366574219e+13},
        {1.559331277284637146e+11,1.149561165512746289e+13,-2.090364987720534766e+13,2.756545234813671484e+13,-3.019842706768454297e+13,2.770540350745921875e+13,-2.176459105054109766e+13,1.225192226734455664e+13,-6.997186266261053467e+10,-1.204215544655462500e+13,2.185911061233155469e+13,-2.703374156063121875e+13,2.938449759477989062e+13,-2.673242217819759375e+13,2.100472118451342969e+13,-1.179704810393715625e+13},
        {1.323100685038486250e+15,-1.360341993819534250e+15,1.376796562643399500e+15,-1.423658727584073000e+15,1.440792882199037250e+15,-1.428285215714143250e+15,1.423892139812061250e+15,-1.434509674775447500e+15,1.435978372553446750e+15,-1.436521328263131500e+15,1.443796679699884000e+15,-1.397169795814290500e+15,1.397732743069400000e+15,-1.367627397169950500e+15,1.367615245144060750e+15,-1.367057611243086750e+15}};
        const double gain13[16] = {1.000000000000000000e0, 9.864749932865497506e-1, 1.008554910983899822e0, 9.909215189952645941e-1, 1.016551297946559407e0, 1.004441025921927633e0, 1.023606422028664920e0, 1.007893430992459427e0,
        1.019020740443200568e0, 1.002733038130677601e0, 1.021101905384611097e0, 1.002615873152301385e0, 1.045142626872535008e0, 9.986874348446839189e-1, 1.022369446580688512e0, 9.918986411869963327e-1};
        const double gain15[16] = {1.000000000000000000e0, 9.927984352486749486e-1, 9.886965960947093901e-1, 9.844666881364172450e-1, 9.811551707487992102e-1, 9.783521372729635512e-1, 9.896903610596464729e-1, 9.869302062198921366e-1,
        9.906593504726489696e-1, 9.859330582601781856e-1, 9.980913383872332956e-1, 9.882277496998294053e-1, 1.000800169181874200e0, 9.924633771768432977e-1, 9.944628454760420233e-1, 1.001591146042800640e0};
        const double align13[2] = {0.7799954861648153, -0.17114181002093012};
        const double align15[2] = {-0.25069580086359916, 1.0363901454673603};
};

int samplelength;
double mean(double *set) {
    double res = 0.0;
    for (int i = 0; i < 52; i++) {
        res += set[i];
    }
    res /= 52;
    return res;
}
double variance(double *set) {
    double res = 0.0;
    double dev[52];
    for (int i = 0; i < 52; i++) {
        dev[i] = pow(set[i]-mean(set), 2);
    }
    res = mean(dev);
    return res;
}
double stadev(double *set) {
    double res = sqrt(variance(set));
    return res;
}
double covar(double *setA, double *setB) {
    double res = 0.0;
    double cov[52];
    for (int i = 0; i < 52; i++) {
        cov[i] = (setA[i] - mean(setA)) * (setB[i] - mean(setB));
    }
    res = mean(cov);
    return res;
}
double corel(double *setA, double *setB) {
    double res = covar(setA, setB)/(stadev(setA)*stadev(setB));
    return res;
}

int ** read_data(char fName[], int address, char mode) {
/*
* COMMANTARY FOR WAVE DATA READING MODE:
*
* 16PickUp Beam Position Monitor(16PU) have a unique readout system, which contain bunch information.
* This function is for reading raw wave binary data from 16pu, and due to the limitation of ADC, the data only have just about 140 turn of bunch.
*
* The structer of wave binary data :
*   Header 19bytes -> wave_0000(year)_00(month)_00(day)_00(hour)_00(min)_00(sec)_address00(13 or 15).dat
*   Wave data [header 4bytes -> wave
*              channel number 2bytes -> 0x00 0xff ~ 0x0f 0xff; no meaning in 0xff
*              wave data 65528*2bytes -> 2bytes is one sample
*              footer 4bytes -> data] 131066btyes -> for 1 channel, and there are 16 channel datas.
* ISSUES
* 1) Channel index of 16PU at #13 of Main Ring(MR): Channel 10 and Channel 11 are swapped.
* 2) ADC of 16PU outputs 14 bits for 1 sample. To make easy to calculate, it has meaningless 0b00 bit so that ADC outputs 16 bits(= 2bytes).
* 3) Can this 140turns limitation be fixed with DMA(Direct Memory Access)?
*
* COMMANTARY FOR PROCESS DATA READING MODE:
*
* 16PU also have a process mode, which contain calculated voltage intensity(or moment).
* With this mode, one bunch have just one sample, ADC can conatin lots of turn. Basicly, process binary data have about 500 turn of bunch.
*
* the structer of process binary data :
*   Header 84bytes -> 0000(year)_00(month)_00(day)_00(hour)_00(min)_00(sec)TRANSVERSE MOMENT measured with sixteen-pu-monitor at address 15
*   Gain setting 2bytes*16 -> original gain * 2^15
*   Data length 3bytes -> unsigned int
*   Process data [data number 2bytes
*                 moment data 3bytes*16 -> 3bytes float]*data length
*   Footer 50bytes -> DATA processed with the 16-pu-monitor circuit
* ISSUES
* 1) Becuase of ADC's 14bits system, processed data is very unique float structer. Unlike float32 or float64, this mode outputs float24,
*    which is consisted of 1bit sign, 5bits exponent, and 18bits fraction.
* 2) At first, it was designed to calculate moment data from bunch, but it is set for just outputing beam intensity with Geortzel algorithm.
*/
    ifstream fin;
    fin.open(fName, ios::in);
    monitor_constant gain;
    if (fin.is_open()) {
        cout << "#" << address << " " << mode << " data file opened." << endl;
        int** voltage = new int*[16];
        if (mode == 0) {
            for (int j = 0; j < 16; j++) {
                voltage[j] = new int[65528];
                fin.seekg(25 + j*131066);
                char buffer[131056];
                fin.read(buffer, 131056);
                for (int i = 0; i < 65528; i++) {
                    if (address == 13) {
                        voltage[j][i] = (((unsigned char) buffer[0 + i*2])*64 + ((unsigned char) buffer[1 + i*2])/4);
                    } else if (address == 15) {
                        voltage[j][i] = (((unsigned char) buffer[0 + i*2])*64 + ((unsigned char) buffer[1 + i*2])/4);
                    }
                }
            }
            for (int i = 0; i < 16; i++) {
                for (int j = 0; j < 65528; j++) {
                    if (address == 13) {
                        voltage[i][j] *= gain.gain13[i];
                    } else if (address == 15) {
                        voltage[i][j] *= gain.gain15[i];
                    }
                } 
            }
        } else if (mode == 1) {
            fin.seekg(116);
            char length_buf[3];
            fin.read(length_buf, 3);
            unsigned int data_length = ((unsigned char) length_buf[0])*65536 + ((unsigned char) length_buf[1])*256 + ((unsigned char) length_buf[2]);
            for (int i = 0; i < data_length; i++) {
                fin.seekg(121 + 50*i);
                char data_buf[48];
                fin.read(data_buf, 48);
                for (int ch = 0; ch < 16; ch++) {
                    voltage[ch] = new int[data_length];
                    bitset<8> bit_buf(data_buf[0+3*ch]);
                    unsigned int exp = ((bit_buf & bitset<8>(124)) >> 2).to_ulong();
                    unsigned int fra = ((bit_buf & bitset<8>(3)).to_ulong())*65536 + ((unsigned char) data_buf[1+3*ch])*256 + ((unsigned char) data_buf[2+3*ch]);
                    double data_buf;
                    double constant = 1/64e6 * 52/64;
                    if (bit_buf[0] == 0) {
                        data_buf = fra * exp * exp * constant;
                    } else {
                        data_buf = (-1) * fra * exp * exp * constant;
                    }
                    if (address == 13) {
                        voltage[ch][i] = data_buf;
                    } else if (address == 15) {
                        voltage[ch][i] = data_buf;
                    }
                }
            }
            for (int i = 0; i < 16; i++) {
                for (int j = 0; j < data_length; j++) {
                    if (address == 13) {
                        voltage[i][j] *= gain.gain13[i];
                    } else if (address == 15) {
                        voltage[i][j] *= gain.gain15[i];
                    }
                } 
            }
        }
        fin.close();
        if (address == 13) {
            int *voltage_buf = new int[65528];
            voltage_buf = voltage[10];
            voltage[10] = voltage[11];
            voltage[11] = voltage_buf;
        }

        return voltage;
    } else {
        cout << "Error opening file." << endl;
        return 0;
    }
}

int * peak_detection(int * voltage) {
/*
* In analysing beam dynamics, the bunch length is also important parameter.
* This function will find start point of beam signal and the bunch length.
* The idea is simple; i+2_th sample is bigger than i+1_th sample, and i+1_th sample is bigger than i_th sample.
* Also i+2_th sample value is bigger than threshold. So the i will be the starting point of first peak.
*/
    int size = 468;
    int num = 0, mean = 0, lag = 0;
    for (int i = 0; i < 65528; i++) {
        num += voltage[i];
    }
    mean = num / 65528;
    int threshold = mean + 50;
    for (int i = 0; i < size; i++) {
        if ((voltage[i] + voltage[i+1] + voltage[i+2])/3 > threshold) {
            lag += i-4;
            break;
        }
    }
    int * voltage_new = new int[65528-lag];
    samplelength = 65528 - lag;
    for (int i = lag; i < 65528; i++) {
        voltage_new[i-lag] = voltage[i];
    }
    return voltage_new;
}

void momentum_calculator(int ** volatage, int bunch, char address) {
/*
* This method is new way to calculate beam bunch moment.
* The paper, which was written by R. Singh for PhD Thesis, introduce the method using least-squre estimating.
* For beam dipole moment, we can assuming that Δ/Σ can be fitted in straigh line.
*
* The fit approach differs fundamentally from the classical integral approach, since it estimates the slope
* parameter (x/r) from the two-dimensional tuples of derived, coupled quantities (Δ,Σ), rather than calculating
* this parameter from the two independent integral estimates of the one-dimensional electrode signals S_L and S_R.
* The measuerment observable is the correlation between the two signals. The intergral approach requires knowledge
* of the baseline to calculate correct integrals, while the fit approach does not!
*/
    monitor_constant matrix_moment;
    double ** result = new double*[16];
    if (address == 13) {
        for (int i = 0; i < 16; i++) {
            result[i] = new double[samplelength];
            for (int j = 0; j < samplelength; j++) {
                for (int k = 0; k < 16; k++) {
                    result[i][j] += matrix_moment.matrix13[i][k] * volatage[k][j];
                }
            }
        }
    } else if (address ==15) {
        for (int i = 0; i < 16; i++) {
            result[i] = new double[samplelength];
            for (int j = 0; j < samplelength; j++) {
                for (int k = 0; k < 16; k++) {
                    result[i][j] += matrix_moment.matrix15[i][k] * volatage[k][j];
                }
            }
        }
    }

    double *sig, *delx, *dely, *x, *y, *sigx, *sigy, *interceptx, *intercepty;
    sig = new double[samplelength];
    delx = new double[samplelength];
    dely = new double[samplelength];
    x = new double[samplelength/52];
    y = new double[samplelength/52];
    sigx = new double[samplelength/52];
    sigy = new double[samplelength/52];
    interceptx = new double[samplelength/52];
    intercepty = new double[samplelength/52];
    for (int i = 0; i < samplelength; i++) {
        sig[i%52] = result[0][i];
        delx[i%52] = result[1][i];
        dely[i%52] = result[2][i];
        if (i % 52 == 51) {
            double dipoleNormal = covar(sig, delx)/variance(sig);
            double dipoleSkew = covar(sig, dely)/variance(sig);
            double sigmaNormal = variance(delx)/(samplelength*variance(sig));
            double sigmaSkew = variance(dely)/(samplelength*variance(sig));
            x[(i+1)/52 - 1] = dipoleNormal;
            y[(i+1)/52 - 1] = dipoleSkew;
            sigx[(i+1)/52 - 1] = sigmaNormal;
            sigy[(i+1)/52 - 1] = sigmaSkew;
            interceptx[(i+1)/52 - 1] = mean(delx) - x[(i+1)/52 - 1]*mean(sig);
            intercepty[(i+1)/52 - 1] = mean(dely) - y[(i+1)/52 - 1]*mean(sig);
            sig = new double[samplelength];
            delx = new double[samplelength];
            dely = new double[samplelength];
        }
    }

    delete [] result;

    ofstream wFile;
    if (address == 13) {
        wFile.open("moment_13.csv");
    } else if (address == 15) {
        wFile.open("moment_15.csv");
    }

    for (int i = 0; i < samplelength/52; i++) {
        if (i % 9 == bunch) {
            string csv = to_string(x[i]) + "," + to_string(sigx[i]) + "," + to_string(interceptx[i]) + "," + to_string(y[i]) + "," + to_string(sigy[i]) + "," + to_string(intercepty[i]) + "," + "\n";
            wFile << csv.c_str();
        }
    }

    wFile.close();

    return;
}

void gnuplot_load() {
/*
* Just for run gnuplot for plotting dipole oscillation.
*/
    FILE *gp;
    gp = popen("gnuplot", "w");
    fprintf(gp, "load \"dipole.gnuplot\"\n");
    pclose(gp);
}

int main(int argc, char **argv) {
/*
* This is a prototype code for read data from 16PU Beam Profile Monitor.
* As a beginer of cpp programming, there may be some wierd coded part for professional, so I kindly ask for understanding.
* Basicly, after compile, drag and drop the #13 and #15 data file on termianl, and press enter to start analysis.
* It use least-square fitting method for calculating beam moment.
*/
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    char wavemode = 0;
    char processmode = 1;

    int **buf_13, **buf_15;
    buf_13 = read_data(argv[1], 13, wavemode);
    buf_15 = read_data(argv[2], 15, wavemode);
    
    int **voltage_13 = new int*[16];
    int **voltage_15 = new int*[16];
    for (int i = 0; i < 16; i++) {
        voltage_13[i] = peak_detection(buf_13[i]);
        voltage_15[i] = peak_detection(buf_15[i]);
    }

    ofstream writevol;
    writevol.open("vol.csv");

    for (int i = 0; i < samplelength; i++) {
        string vol_csv = to_string(voltage_13[0][i]) + "," + to_string(voltage_13[1][i]) + "," + to_string(voltage_13[2][i]) + "," + to_string(voltage_13[3][i]) + "\n";
        writevol << vol_csv.c_str();
    }

    writevol.close();

    delete [] buf_13;
    delete [] buf_15;

    int bunch;
    cout << " Select bunch number : ";
    cin >> bunch;

    momentum_calculator(voltage_13, bunch, 13);
    momentum_calculator(voltage_15, bunch, 15);

    end = chrono::system_clock::now();
    const double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count()/1000000.0);
    cout << "Analysis Done." << endl;
    cout << "Execution time : " << time << "s" << endl;

    gnuplot_load();

    return 0;
}
