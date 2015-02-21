#include <iostream>
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/core/core.hpp"
#include <string>
#include <vector>

using namespace std;
using namespace cv;

struct piece 
{
    int id;
    int length;
    int width;
    int depth;
    int volume;
};

/*******************************************
* This function takes a picture using a
* webcam. Then stores the picture to the 
* specified file location. Device 
* number is passed as 'cam'. 
*******************************************/
void imageCap(int cam, string fileLocation)
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
    
    //imwrite("/Users/vegas_bballer/Documents/Senior_Project/image.jpg",frame);
    
    if (cam == 0)
    {
        imwrite(fileLocation, frame);
        cout << "Right side written to a file.\n";
    }
    else if (cam == 1)
    {
        imwrite(fileLocation, frame);
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
* This function figures out if there are
* more pieces to be processed in the image
*******************************************/
bool moreObjects()
{
    bool more = false;
    
    return more; 
}

/*******************************************
* This function finds piece length
*******************************************/
int findPieceLength()
{
    int l;
    
    return l;
}

/*******************************************
* This function finds piece width
*******************************************/
int findPieceWidth()
{
    int w;
    
    return w;
}

/*******************************************
* This function finds piece depth
*******************************************/
int findPieceDepth()
{
    int d;
    
    return d;
}

/*******************************************
* This function processes the left and right
* pictures storing measurments in a vector
*******************************************/
vector process(string rightFile, string leftFile)
{
    vector<piece> sample;
    piece zero;
    zero.length = 0;
    zero.width = 0;
    zero.depth = 0;
    zero.volume = 0;
    
    bool finished = false;
    int i = 0;
    
    //loop through all pieces
    while (!finished)
    {
        sample.insert(0,zero);
        sample[0].id = i;
        sample[0].length = findPieceLength();
        sample[0].width = findPieceWidth();
        sample[0].depth = findPieceDepth();
        sample[0].volume = sample[0].length * sample[0].width * sample[0].depth;
        finished = !moreObjects();
        i++;
    }
    return sample;
}

/*******************************************
* This function uploads the measurments to
* an excel file 
*******************************************/
void upload(string uploadFile, vector sample)
{
    
}

/*******************************************
* Main
*******************************************/
int main(int argc, const char * argv[])
{
    string rightFile;
    string leftFile;
    cout << "Right pic file location:";
    cin >> rightFile;
    cout << "Left pic file location:";
    cin leftFile;
    
    int rightCam = 0;
    int leftCam = 1;
    
    imageCap(rightCam, rightFile);
    imageCap(leftCam, leftFile;

    //video(rightCam, leftCam);
    
    string uploadFile;
    cout << "Upload measurments to:"
    cin >> uploadFile;
    upload(uploadFile, process(rightFile, leftFile););
    
    return 0;
}
