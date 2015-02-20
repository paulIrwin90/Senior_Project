#include <iostream>
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/core/core.hpp"

using namespace std;
using namespace cv;

/*******************************************
* This function takes a picture using a
* webcam. Then stores the picture. Device 
* number is passed as 'cam'. More 
* processing can be done in this function.
*******************************************/
void imageCap(int cam)
{
    VideoCapture webcam(cam);
    Mat frame, grayed, newFrame;
    Size size(960, 540);
    
    webcam.read(frame);
    resize(frame, frame, size);
    
    //cvtColor(frame, grayed, CV_RGB2GRAY);
    //blur(grayed, newFrame, Size(6, 6));
    //Canny(newFrame, newFrame, 30, 10);
    
    imshow("Picture", frame);
    waitKey(0);             // Press a key to proceed
    
    destroyAllWindows();    // clean up screen
    
    if (cam == 1)
    {
        imwrite("/Users/vegas_bballer/Documents/rightSide.jpg", newFrame);
        cout << "Right side written to a file.\n";
    }
    else if (cam == 2)
    {
        imwrite("/Users/vegas_bballer/Documents/leftSide.jpg", newFrame);
        cout << "Left side written to a file.\n";
    }
    else
        cout << "Error: Device number invalid.\n";
}

/*******************************************
* This function captures a video using both
* cameras simultaneously. Use 'esc' to end
* video stream.
*******************************************/
void video(int cam1, int cam2)
{
    VideoCapture webcam1(cam1), webcam2(cam2);
    Mat video1, video2;
    
    webcam1.set(CV_CAP_PROP_FRAME_HEIGHT, 540);
    webcam1.set(CV_CAP_PROP_FRAME_WIDTH, 540);
    
    webcam2.set(CV_CAP_PROP_FRAME_HEIGHT, 540);
    webcam2.set(CV_CAP_PROP_FRAME_WIDTH, 540);
    
    while(true)
    {
        webcam1.read(video1);
        imshow("Video 1", video1);
        webcam2.read(video2);
        imshow("Video 2", video2);
        
        if (waitKey(30) == 27)          // Press esc to proceed
        {
            destroyAllWindows();        // clear up screen
            break;
        }
    }
}

/*******************************************
* Main
*******************************************/
int main(int argc, const char * argv[])
{
    int rightCam = 0;
    int leftCam = 1;
    
    imageCap(rightCam);
    imageCap(leftCam);

    video(rightCam, leftCam);
    
    return 0;
}
