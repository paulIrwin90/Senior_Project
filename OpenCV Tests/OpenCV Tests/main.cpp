#include <iostream>
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/core/core.hpp"

using namespace std;
using namespace cv;

int imageCap(int cam)
{
    VideoCapture webcam(cam);
    Mat frame;
    
    webcam.read(frame);
    
    imshow("Picture", frame);
    waitKey(0);
    //imwrite("/Users/vegas_bballer/Documents/new.jpg", frame);
    
    return 0;
}

int main(int argc, const char * argv[])
{
    imageCap(0);
    imageCap(1);
    return 0;
}
