/*************************************************************
* INL BIO MATERIAL VOLUME CAMERA SENSOR
*
* This code controls the operation of the INL BIO MATERIAL
* VOLUME CAMERA SENSOR. It was written by
* -Quinn Stratton
* -Charles Hasler
* -Paul Irwin
* -Cody Marshall
* -David Fluckiger
*
* The code operates two cameras mounted side-by-side over a
* conveyor belt. The conveyor belt holds wood chips that need
* to be measured by the system. Length, width, area, depth,
* and volume are the main outputs being stored.
*
* Data is output to images an Excel spreadsheets. 
*************************************************************/

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/contrib/contrib.hpp"

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
//#include <windows.h>
#include <ctime>
//#include <direct.h>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <ctype.h>

using namespace cv;
using namespace std;

/*************************************************************
* Struct of Each Measured Object
*************************************************************/
struct piece {
    int pictureID;
    int objectID;
    double length;
    double width;
    double perimeter;
    double area;
    double depth;
    double volume;
    double posX;
    double posY;
    
};


/*************************************************************
* Global Variables
*************************************************************/
const int RIGHT_CAM = 0;
const int LEFT_CAM = 1;
const int THRESHOLD = 160;
const int CANNY_THRESH = 250;
const int MAX_THRESH = 255;
RNG rng(12345);
int CURRENT_ID = 0;
int PICTURE_ID = 0;
string dir;
string imgDir;

/*************************************************************
* Function Declarations
*************************************************************/

//OUTPUT calibration();

int calibration(string intrinsicFile, string extrinsicFile);
void stereoCalibCapture(Size boardsize);
bool readStringList( const string& filename, vector<string>& l );
bool stereoCalib(const vector<string>& imagelist, Size boardSize, string intrinsicFile, string extrinsicFile);
int stereoMatch(Mat imgL, Mat imgR, string intrinsicFile, string extrinsicFile);


void take_picture(Mat &img, int cam);
void detect_and_measure(Mat &imgL, Mat &imgR, vector<piece> &piecesVec, string path);
void export_data(vector<piece> &pz, string fileName, bool append);
bool dirExists(const string&);
string getFileName();
string intToString(int integer, int width, char fill);

/*************************************************************
* Main Function
*************************************************************/
int main(int argc, char** argv)
{
    cout << "Initializing INL BIO MATERIAL VOLUME CAMERA SENSOR" << endl;
    cout << "Please enter a save location directory: ";
    string save_location;
    cin >> save_location;
    string fileName;
    

    // Checks to see if this directory exists, if not it asks to create it
    if (dirExists(save_location)) {
        fileName = getFileName();
        cout << "These measurements will be saved in " << fileName << endl;
    }
    else {
        cout << "That directory doesn't exist. Would you like it created? (y/n): ";
        char ans;
        cin >> ans;
        if ((ans == 'y') || (ans == 'Y')) {
            if (!(_mkdir(save_location.c_str()))) {
                cout << "Error in creating directory" << endl;
                exit(EXIT_FAILURE);
            }
        }
        else if ((ans == 'n') || (ans == 'N')) {
            cout << "Exiting..." << endl;
            exit(EXIT_SUCCESS);
        }
        else {
            cout << "Option not recognized." << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    // Tells the user what the excel file will be called
    cout << "Measurements will be published to " << fileName << endl;
    
    // Declares needed variables
    vector<piece> pieces;
    
    // Starts Calibration
    string intrinsicFile = dir + "intrinsics.yml";
    string extrinsicFile = dir + "extrinsics.yml";
    calibration(intrinsicFile, extrinsicFile);
    
    // Takes the two pictures and saves them to a temporary file
    Mat left_ori_img;
    Mat right_ori_img;

    bool done = false;
    do
    {
        take_picture(left_ori_img, LEFT_CAM);
        take_picture(right_ori_img, RIGHT_CAM);
    
        // Gets the zero-to-threshold images
        Mat left_zero = left_ori_img;
        Mat right_zero = right_ori_img;
    
        // Performs stereo matching for disparity map
        stereoMatch(left_ori_img, right_ori_img, intrinsicFile, extrinsicFile);
    
        // Detects objects and records measurements into a struct
        detect_and_measure(left_ori_img, right_ori_img, pieces, save_location);
        
    } while (!done);

    // Exports the data to an Excel file
    export_data(pieces, (save_location + "/" + fileName), false);
    
    

    waitKey(0);
    
    return(0);
}


//===========================================
//
//
//===========================================
int calibration(string intrinsicFile, string extrinsicFile)
{
    bool calibrated = false;
    // Corners of Chessboard width y height
    Size boardSize = Size(9, 6);
    
    // Capture the images and save to imgDir
    stereoCalibCapture(boardSize);
    
    // Read file with pictures
    string imagelistfn = dir + "stereo_calib.xml";
    vector<string> imagelist;
    
    imagelist.resize(0);
    FileStorage fs(imagelistfn, FileStorage::READ);
    if( !fs.isOpened() ) {
        cout << "Error opening file " << imagelistfn << endl;
        return -1;
    }
    
    FileNode n = fs.getFirstTopLevelNode();
    if( n.type() != FileNode::SEQ ) {
        cout << "File type is incorrect.\n";
        return -1;
    }
    
    FileNodeIterator it = n.begin(), it_end = n.end();
    for( ; it != it_end; ++it ) {
        imagelist.push_back((string)*it);
    }
    
    if (imagelist.empty()) {
        cout << "The string list is empty." << endl;
        return -1;
    }
    else if( imagelist.size() % 2 != 0 ) {
        cout << "Error: the image list contains odd (non-even) number of elements\n";
        return -1;
    }
    else // Start Calibration
        calibrated = stereoCalib(imagelist, boardSize, intrinsicFile, extrinsicFile);
    
    if (!calibrated) {
        cout << "There was an error durning calibration.\n";
        return -1;
    }
    else
        cout << "Stereo Calibration complete.\n";
    
    return 0;
}


//===========================================
//
//
//===========================================
void stereoCalibCapture(Size boardsize)
{
    VideoCapture capR(RIGHT_CAM), capL(LEFT_CAM);
    Mat right, left, grayR, grayL;
    Mat cornersR, cornersL;
    
    int k, count = 0, numpairs = 6;
    bool foundR = false, foundL = false;
    
    // Example file location:
    // C:/opencv/images/
    string fileR, fileL;
    string path = dir + imgDir + "image";
    string r = "Right";
    string l = "Left";
    string ext = ".jpg";
    
    // Set height and width of images
    //capR.set(CV_CAP_PROP_FRAME_HEIGHT, 700);
    //capR.set(CV_CAP_PROP_FRAME_WIDTH, 800);
    
    //capL.set(CV_CAP_PROP_FRAME_HEIGHT, 700);
    //capL.set(CV_CAP_PROP_FRAME_WIDTH, 800);
    
    cout << "Starting calibration image capture.\n"
         << "How many pairs of images would you like to capture? ";
    cin >> numpairs;
    cout << "To save image press 'return'.  -  To reset image press 'esc'.";
        
    while(count < numpairs)
    {
        capR >> right;
        capL >> left;
        
        // Change image size
        //resize(right, right, Size(800, 700));
        //resize(left, left, Size(800, 700));
        
        cvtColor(right, grayR, CV_RGB2GRAY);
        cvtColor(left, grayL, CV_RGB2GRAY);
        
        // Test pictures to see if the chessboard corners can be found
        foundR = findChessboardCorners(right, boardsize, cornersR, CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_NORMALIZE_IMAGE);
        foundL = findChessboardCorners(left, boardsize, cornersL, CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_NORMALIZE_IMAGE);
        
        // Show viewer the corners found
        if (foundR)
        {
            cornerSubPix(grayR, cornersR, Size(11, 11), Size(-1, -1), TermCriteria(CV_TERMCRIT_EPS | CV_TERMCRIT_ITER, 30, 0.1));
            drawChessboardCorners(grayR, boardsize, cornersR, foundR);
        }
        if (foundL)
        {
            cornerSubPix(grayL, cornersL, Size(11, 11), Size(-1, -1), TermCriteria(CV_TERMCRIT_EPS | CV_TERMCRIT_ITER, 30, 0.1));
            drawChessboardCorners(grayL, boardsize, cornersL, foundL);
        }
        
        // Show images
        imshow("Calibrate Right", grayR);
        imshow("Calibrate Left", grayL);
        
        k = waitKey(10);
        
        if (k == 13 )//&& foundR && foundL) // return
        {
            cout << "Found corners. Saving images.\n";
            fileR = path + r + to_string(count + 1) + ext;
            fileL = path + l + to_string(count + 1) + ext;
            
            imwrite(fileR, right);
            imwrite(fileL, left);
            count++;
        }
        
        else if (k == 27) // esc
        {
            cout  << "Retrying image capture.\n";
            destroyAllWindows();
        }
    }
    destroyAllWindows();
}


//===========================================
//
//
//===========================================
bool stereoCalib(const vector<string>& imagelist, Size boardSize, string intrinsicFile, string extrinsicFile)
{
    bool displayCorners = true;
    const int maxScale = 2;
    const float squareSize = 1.f;  // Set this to your actual square size
    // ARRAY AND VECTOR STORAGE:
    
    vector<vector<Point2f> > imagePoints[2];
    vector<vector<Point3f> > objectPoints;
    Size imageSize;
    
    int i, j, k, nimages = (int)imagelist.size()/2;
    
    imagePoints[0].resize(nimages);
    imagePoints[1].resize(nimages);
    vector<string> goodImageList;
    
    for( i = j = 0; i <= nimages; i++ )
    {
        for( k = 0; k < 2; k++ )
        {
            const string& filename = imagelist[i*2+k];
            Mat img = imread(filename, 0);
            if(img.empty())
                break;
            if( imageSize == Size() )
                imageSize = img.size();
            else if( img.size() != imageSize )
            {
                cout << "The image " << filename << " has the size different from the first image size. Skipping the pair\n";
                break;
            }
            bool found = false;
            vector<Point2f>& corners = imagePoints[k][j];
            for( int scale = 1; scale <= maxScale; scale++ )
            {
                Mat timg;
                if( scale == 1 )
                    timg = img;
                else
                    resize(img, timg, Size(), scale, scale);
                    found = findChessboardCorners(timg, boardSize, corners,
                                              CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_NORMALIZE_IMAGE);
                if( found ) {
                    if( scale > 1 ) {
                        Mat cornersMat(corners);
                        cornersMat *= 1./scale;
                    }
                    break;
                }
            }
            if( displayCorners ) {
                cout << "Showing " << filename << endl;
                Mat cimg, cimg1;
                cvtColor(img, cimg, COLOR_GRAY2BGR);
                drawChessboardCorners(cimg, boardSize, corners, found);
                double sf = 640./MAX(img.rows, img.cols);
                resize(cimg, cimg1, Size(), sf, sf);
                imshow("Chessboard Corners", cimg1);
                waitKey(0);
            }
            else
                putchar('.');
            if( !found )
                break;
            cornerSubPix(img, corners, Size(11,11), Size(-1,-1),
                         TermCriteria(CV_TERMCRIT_ITER+CV_TERMCRIT_EPS, 30, 0.01));
        }
        if( k == 2 ) {
            goodImageList.push_back(imagelist[i*2]);
            goodImageList.push_back(imagelist[i*2+1]);
            j++;
        }
    }
    cout << j << " pairs have been successfully detected.\n";
    nimages = j;
    if( nimages < 2 ) {
        cout << "Error: too little pairs to run the calibration\n";
        return false;
    }
    
    imagePoints[0].resize(nimages);
    imagePoints[1].resize(nimages);
    objectPoints.resize(nimages);
    
    for( i = 0; i < nimages; i++ )
    {
        for( j = 0; j < boardSize.height; j++ )
            for( k = 0; k < boardSize.width; k++ )
                objectPoints[i].push_back(Point3f(j*squareSize, k*squareSize, 0));
    }
    
    cout << "Running stereo calibration ...\n";
    
    Mat cameraMatrix[2], distCoeffs[2];
    cameraMatrix[0] = Mat::eye(3, 3, CV_64F);
    cameraMatrix[1] = Mat::eye(3, 3, CV_64F);
    Mat R, T, E, F;
    
    double rms = stereoCalibrate(objectPoints, imagePoints[0], imagePoints[1],
                                 cameraMatrix[0], distCoeffs[0],
                                 cameraMatrix[1], distCoeffs[1],
                                 imageSize, R, T, E, F,
                                 TermCriteria(CV_TERMCRIT_ITER+CV_TERMCRIT_EPS, 100, 1e-5),
                                 CV_CALIB_FIX_ASPECT_RATIO +
                                 CV_CALIB_ZERO_TANGENT_DIST +
                                 CV_CALIB_SAME_FOCAL_LENGTH +
                                 CV_CALIB_RATIONAL_MODEL +
                                 CV_CALIB_FIX_K3 + CV_CALIB_FIX_K4 + CV_CALIB_FIX_K5);
    cout << "done with RMS error=" << rms << endl;
    
    // CALIBRATION QUALITY CHECK
    // because the output fundamental matrix implicitly
    // includes all the output information,
    // we can check the quality of calibration using the
    // epipolar geometry constraint: m2^t*F*m1=0
    double err = 0;
    int npoints = 0;
    vector<Vec3f> lines[2];
    for( i = 0; i < nimages; i++ )
    {
        int npt = (int)imagePoints[0][i].size();
        Mat imgpt[2];
        for( k = 0; k < 2; k++ )
        {
            imgpt[k] = Mat(imagePoints[k][i]);
            undistortPoints(imgpt[k], imgpt[k], cameraMatrix[k], distCoeffs[k], Mat(), cameraMatrix[k]);
            computeCorrespondEpilines(imgpt[k], k+1, F, lines[k]);
        }
        for( j = 0; j < npt; j++ )
        {
            double errij = fabs(imagePoints[0][i][j].x*lines[1][j][0] +
                                imagePoints[0][i][j].y*lines[1][j][1] + lines[1][j][2]) +
            fabs(imagePoints[1][i][j].x*lines[0][j][0] +
                 imagePoints[1][i][j].y*lines[0][j][1] + lines[0][j][2]);
            err += errij;
        }
        npoints += npt;
    }
    cout << "average reprojection err = " <<  err/npoints << endl;
    
    // save intrinsic parameters
    FileStorage fs(intrinsicFile, CV_STORAGE_WRITE);
    if( fs.isOpened() ) {
        fs << "M1" << cameraMatrix[0] << "D1" << distCoeffs[0] <<
        "M2" << cameraMatrix[1] << "D2" << distCoeffs[1];
        fs.release();
    }
    else
        cout << "Error: can not save the intrinsic parameters\n";
    
    Mat R1, R2, P1, P2, Q;
    Rect validRoi[2];
    
    stereoRectify(cameraMatrix[0], distCoeffs[0],
                  cameraMatrix[1], distCoeffs[1],
                  imageSize, R, T, R1, R2, P1, P2, Q,
                  CALIB_ZERO_DISPARITY, 1, imageSize, &validRoi[0], &validRoi[1]);
    
    fs.open(extrinsicFile, CV_STORAGE_WRITE);
    if( fs.isOpened() )
    {
        fs << "R" << R << "T" << T << "R1" << R1 << "R2" << R2 << "P1" << P1 << "P2" << P2 << "Q" << Q;
        fs.release();
    }
    else
        cout << "Error: can not save the extrinsic parameters\n";
    
    // Compute and display rectification
    Mat rmap[2][2];
    
    // Calibration using Hartley's method
    // use intrinsic parameters of each camera, but
    // compute the rectification transformation directly
    // from the fundamental matrix

    vector<Point2f> allimgpt[2];
    for( k = 0; k < 2; k++ ) {
        for( i = 0; i < nimages; i++ )
            std::copy(imagePoints[k][i].begin(), imagePoints[k][i].end(), back_inserter(allimgpt[k]));
    }
    
    F = findFundamentalMat(Mat(allimgpt[0]), Mat(allimgpt[1]), FM_8POINT, 0, 0);
    Mat H1, H2;
    stereoRectifyUncalibrated(Mat(allimgpt[0]), Mat(allimgpt[1]), F, imageSize, H1, H2, 3);
    
    R1 = cameraMatrix[0].inv()*H1*cameraMatrix[0];
    R2 = cameraMatrix[1].inv()*H2*cameraMatrix[1];
    P1 = cameraMatrix[0];
    P2 = cameraMatrix[1];

    //Precompute maps for cv::remap()
    initUndistortRectifyMap(cameraMatrix[0], distCoeffs[0], R1, P1, imageSize, CV_16SC2, rmap[0][0], rmap[0][1]);
    initUndistortRectifyMap(cameraMatrix[1], distCoeffs[1], R2, P2, imageSize, CV_16SC2, rmap[1][0], rmap[1][1]);
    
    Mat canvas;
    double sf;
    int w, h;
    
    sf = 300./MAX(imageSize.width, imageSize.height);
    w = cvRound(imageSize.width*sf);
    h = cvRound(imageSize.height*sf);
    canvas.create(h*2, w, CV_8UC3);

    
    for( i = 0; i < nimages; i++ )
    {
        for( k = 0; k < 2; k++ )
        {
            Mat img = imread(goodImageList[i*2+k], 0), rimg, cimg;
            remap(img, rimg, rmap[k][0], rmap[k][1], CV_INTER_LINEAR);
            cvtColor(rimg, cimg, COLOR_GRAY2BGR);

            Mat canvasPart = canvas(Rect(w*k, 0, w, h));
            resize(cimg, canvasPart, canvasPart.size(), 0, 0, CV_INTER_AREA);
        }
        
        // Epipolar lines
        for( j = 0; j < canvas.rows; j += 16 )
            line(canvas, Point(0, j), Point(canvas.cols, j), Scalar(0, 255, 0), 1, 8);
        
        
        imshow("Rectified", canvas);
        string file = dir + "rectified" + to_string(i + 1) + ".jpg";
        imwrite(file, canvas);
        
        waitKey(0);
        
    }
}


//===========================================
//
//
//===========================================
int stereoMatch(Mat imgL, Mat imgR, string intrinsicFile, string extrinsicFile)
{
    string disparity_filename = dir + "disparity.jpg";
    string point_cloud_filename = dir + "pointcloud.jpg";
    
    // These are the possible algorithms
    enum
    {
        STEREO_BM   = 0,
        STEREO_SGBM = 1,
        STEREO_HH   = 2,
        STEREO_VAR  = 3
    };
    
    int alg = STEREO_SGBM;
    int SADWindowSize = 19; // this is the block size, must be positive and odd...
    int numberOfDisparities = 11 * 16; // must be divisible by 16
    bool no_display = false;
    float scale = 1.0; // must be positive floating point
    
    StereoBM bm;
    StereoSGBM sgbm;
    StereoVar var;
    
    int color_mode = alg == STEREO_BM ? 0 : -1; // if algorithm is StereoBM use 0 otherwise -1
    if (alg != STEREO_SGBM) {
        cvtColor(imgL, imgL, color_mode);
        cvtColor(imgR, imgR, color_mode);
    }
    
    Size img_size = imgL.size();
    Rect roi1, roi2;
    Mat Q;
    
    // Get intrinsic and extrinsic parameters from files
    if( intrinsicFile != "" )
    {
        // reading intrinsic parameters
        FileStorage fs(intrinsicFile, CV_STORAGE_READ);
        if(!fs.isOpened())
        {
            cout << "Failed to open file " << intrinsicFile << endl;
            return -1;
        }
        
        Mat M1, D1, M2, D2;
        fs["M1"] >> M1;
        fs["D1"] >> D1;
        fs["M2"] >> M2;
        fs["D2"] >> D2;
        
        M1 *= scale;
        M2 *= scale;
        
        // reading extrinsic parameters
        fs.open(extrinsicFile, CV_STORAGE_READ);
        if(!fs.isOpened())
        {
            cout << "Failed to open file " << extrinsicFile << endl;
            return -1;
        }
        
        Mat R, T, R1, P1, R2, P2;
        fs["R"] >> R;
        fs["T"] >> T;
        
        stereoRectify( M1, D1, M2, D2, img_size, R, T, R1, R2, P1, P2, Q, CALIB_ZERO_DISPARITY, -1, img_size, &roi1, &roi2 );
        
        Mat map11, map12, map21, map22;
        initUndistortRectifyMap(M1, D1, R1, P1, img_size, CV_16SC2, map11, map12);
        initUndistortRectifyMap(M2, D2, R2, P2, img_size, CV_16SC2, map21, map22);
        
        Mat imgLr, imgRr;
        remap(imgL, imgLr, map11, map12, INTER_LINEAR);
        remap(imgR, imgRr, map21, map22, INTER_LINEAR);
        
        imgL = imgLr;
        imgR = imgRr;
    }
    
    numberOfDisparities = numberOfDisparities > 0 ? numberOfDisparities : ((img_size.width/8) + 15) & -16;
    
    // Depending on which algorithm used
    bm.state->roi1 = roi1;
    bm.state->roi2 = roi2;
    bm.state->preFilterCap = 31;
    bm.state->SADWindowSize = SADWindowSize > 0 ? SADWindowSize : 15;
    bm.state->minDisparity = 0;
    bm.state->numberOfDisparities = numberOfDisparities;
    bm.state->textureThreshold = 10;
    bm.state->uniquenessRatio = 15;
    bm.state->speckleWindowSize = 0;
    bm.state->speckleRange = 1;
    bm.state->disp12MaxDiff = 1;
    
    sgbm.preFilterCap = 63;
    sgbm.SADWindowSize = SADWindowSize > 0 ? SADWindowSize : 3;
    int cn = imgL.channels();
    sgbm.P1 = 8 * cn * sgbm.SADWindowSize * sgbm.SADWindowSize;
    sgbm.P2 = 32 * cn * sgbm.SADWindowSize * sgbm.SADWindowSize;
    sgbm.minDisparity = 0;
    sgbm.numberOfDisparities = numberOfDisparities;
    sgbm.uniquenessRatio = 30;
    sgbm.speckleWindowSize = bm.state->speckleWindowSize;
    sgbm.speckleRange = bm.state->speckleRange;
    sgbm.disp12MaxDiff = 1;
    sgbm.fullDP = alg == STEREO_SGBM;
    
    // Recommended from other source... not tested yet
    //sgbm.SADWindowSize = 5;
    //sgbm.numberOfDisparities = 192;
    //sgbm.preFilterCap = 4;
    //sgbm.minDisparity = -64;
    //sgbm.uniquenessRatio = 1;
    //sgbm.speckleWindowSize = 150;
    //sgbm.speckleRange = 2;
    //sgbm.disp12MaxDiff = 10;
    //sgbm.fullDP = false;
    //sgbm.P1 = 600;
    //sgbm.P2 = 2400;
    
    var.levels = 3;                                 // ignored with USE_AUTO_PARAMS
    var.pyrScale = 0.5;                             // ignored with USE_AUTO_PARAMS
    var.nIt = 25;
    var.minDisp = -numberOfDisparities;
    var.maxDisp = 0;
    var.poly_n = 3;
    var.poly_sigma = 0.0;
    var.fi = 15.0f;
    var.lambda = 0.03f;
    var.penalization = var.PENALIZATION_TICHONOV;   // ignored with USE_AUTO_PARAMS
    var.cycle = var.CYCLE_V;                        // ignored with USE_AUTO_PARAMS
    var.flags = var.USE_SMART_ID | var.USE_AUTO_PARAMS | var.USE_INITIAL_DISPARITY | var.USE_MEDIAN_FILTERING ;
    
    Mat disp, disp8;
    //Mat img1p, img2p, dispp;
    //copyMakeBorder(img1, img1p, 0, 0, numberOfDisparities, 0, IPL_BORDER_REPLICATE);
    //copyMakeBorder(img2, img2p, 0, 0, numberOfDisparities, 0, IPL_BORDER_REPLICATE);
    
    if( alg == STEREO_BM )
        bm(imgL, imgR, disp);
    else if( alg == STEREO_VAR ) {
        var(imgL, imgR, disp);
    }
    else if( alg == STEREO_SGBM || alg == STEREO_HH )
        sgbm(imgL, imgR, disp);
    
    
    //disp = dispp.colRange(numberOfDisparities, img1p.cols);
    if( alg != STEREO_VAR )
        disp.convertTo(disp8, CV_8U, 255/(numberOfDisparities*16.));
    else
        disp.convertTo(disp8, CV_8U);
    
    if( !no_display )
    {
        //namedWindow("left", 1);
        imshow("left", imgL);
        //namedWindow("right", 1);
        imshow("right", imgR);
        //namedWindow("disparity", 0);
        imshow("disparity", disp8);
        
        cout << "Press any key to continue...\n";
        waitKey(0);
    }
    
    if(disparity_filename != "")
        imwrite(disparity_filename, disp8);
    
    // This does not work yet
    //if(point_cloud_filename != "")
    //{
        //printf("storing the point cloud...");
        //fflush(stdout);
        //Mat xyz;
        //reprojectImageTo3D(disp, xyz, Q, true);
        //saveXYZ(point_cloud_filename, xyz);
        //printf("\n");
    //}
    
    return 0;
}



//===========================================
// Function: dirExists
// Description: Checks to see if a given
// directory exists.
// Parameters: string& dirName_in - The
// path to the directory
// Returns: bool - true if exists
//===========================================
bool dirExists(const string& dirName_in)
{
    DWORD ftyp = GetFileAttributesA(dirName_in.c_str());
    if (ftyp == INVALID_FILE_ATTRIBUTES)
        return false; //something is wrong with your path!

    if (ftyp & FILE_ATTRIBUTE_DIRECTORY)
        return true; // this is a directory!
    
    return false; // this is not a directory!
    
}


//===========================================
// Function: createFileName
// Description: Creates file name base on
// date and time
// Parameters: NONE
// Returns: string - Excel file name
//===========================================
string getFileName()
{
    time_t t = time(0); // get time now
    struct tm * now = localtime(&t);
    string mon = intToString(now->tm_mon + 1, 2, '0');
    string day = intToString(now->tm_mday, 2, '0');
    string hour = intToString(now->tm_hour, 2, '0');
    string min = intToString(now->tm_min, 2, '0');
    string sec = intToString(now->tm_sec, 2, '0');
    
    return (to_string(now->tm_year + 1900) + "_" + mon + "_" + day + "_" + hour + "_" + min + "_" + sec + ".csv");
    
}


//===========================================
// Function: intToString
// Description: Converts an integer to a
// string with setw and setfill
// Parameters: int integer - int to be
// converted
// int width - set width value
// cahr fill - set fill char
// Returns: string - converted int
//===========================================
string intToString(int integer, int width, char fill)
{
    ostringstream ss;
    ss << setw(width) << setfill(fill) << integer;
    return ss.str();
}


//===========================================
// Function: takePicture
// Description: Takes a picture and stores it
// in the img Mat
// Parameters: Mat img - Mat to hold capture
// int cam - 0 for camera 1
// 1 for camera 2
// Returns: none
//===========================================
void take_picture(Mat &img, int cam)
{
    VideoCapture webcam(cam);
    if (!webcam.isOpened())	{
        cout << "Error: Unable to communicate with camera" << endl;
    }
    else {
        webcam >> img;
    }
    
}


//===========================================
// Function: detect_and_measure
// Description: Looks at an image, identifies
// all the pieces, creates a
// struct for each one, and
// records measurements for each
// Parameters: Mat imgL - orig. left image
// Mat imgR - orig. right image
// vector<piece> - vector or
// pieces
// Returns: none
//===========================================
void detect_and_measure(Mat &imgL, Mat &imgR, vector<piece> &piecesVec, string path)
{
    // Creates Window
    char* source_window = "Source";
    namedWindow(source_window, CV_WINDOW_AUTOSIZE);
    imshow(source_window, imgL);

    // Convert image to gray and blur it
    Mat img_mod;
    Mat zero_modL;
    Mat zero_modR;

    cvtColor(imgL, img_mod, CV_BGR2GRAY);
    blur(img_mod, img_mod, Size(3, 3));
    
    // Gets a binary threshold (Black and White Image)
    threshold(img_mod, img_mod, THRESHOLD, MAX_THRESH, THRESH_BINARY);
    
    // Gets a zero-to-threshold (Background turned to black)
    blur(imgL, zero_modL, Size(3, 3));
    blur(imgR, zero_modR, Size(3, 3));
    threshold(imgL, zero_modL, THRESHOLD, MAX_THRESH, THRESH_TOZERO_INV);
    threshold(imgR, zero_modR, THRESHOLD, MAX_THRESH, THRESH_TOZERO_INV);
    
    // Saves zero-to-threshold left and right images
    string full_pathL = path + "/Picture_Left_" + to_string(PICTURE_ID) + ".jpg";
    imwrite(full_pathL, zero_modL);
    
    string full_pathR = path + "/Picture_Right_" + to_string(PICTURE_ID) + ".jpg";
    imwrite(full_pathR, zero_modR);
    
    // Detects edges
    vector<vector<Point> > contours;
    vector<Vec4i> hierarchy;
    
    Canny(img_mod, img_mod, CANNY_THRESH, 2 * CANNY_THRESH, 5);
    img_mod.convertTo(img_mod, CV_8U);
    
    // Find contours
    findContours(img_mod, contours, hierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
    
    // Find the rotated rectangles and ellipses for each contour
    vector<RotatedRect> minRect(contours.size());
    int largest_contour = 0;

    // Draw contours + rotated rects
    Mat drawing = Mat::zeros(img_mod.size(), CV_8UC3);
    Scalar color;
    for (int i = 0; i < contours.size(); i++)
    {
        minRect[i] = minAreaRect(Mat(contours[i]));
        // Changes the color for each piece
        color = Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
        
        // Draws contours on original image
        drawContours(imgL, contours, i, color, 3, 8, vector<Vec4i>(), 0, Point());
        
        // Draws rectangle around piece
        Point2f rect_points[4];
        minRect[i].points(rect_points);
        for (int j = 0; j < 4; j++) {
            line(imgL, rect_points[j], rect_points[(j + 1) % 4], color, 1, 8);
        }
        
        piece pz;
        pz.length = minRect[i].size.height;
        pz.width = minRect[i].size.width;
        pz.area = contourArea(contours[i], false);
        pz.perimeter = arcLength(contours[i], true);
        pz.objectID = CURRENT_ID++;
        pz.posX = minRect[i].center.x;
        pz.posY = minRect[i].center.y;
        pz.pictureID = PICTURE_ID;
        piecesVec.push_back(pz);
        
        // Prints OBJECT_ID on the image
        putText(imgL, to_string(pz.objectID), minRect[i].center, FONT_HERSHEY_COMPLEX_SMALL, 0.5, CV_RGB(0, 0, 0), 1.0);
    }
    
    // Saves image
    string full_path = path + "/Picture_" + to_string(PICTURE_ID) + ".jpg";
    imwrite(full_path, imgL);
    
    // Increments PICTURE_ID for next image
    PICTURE_ID++;
}


//===========================================
// Function: export_data
// Description: Prints all of the measured
// data in an excel sheet
// Parameters:
// Returns: none
//===========================================
void export_data(vector<piece> &pieces, string fileName, bool append)
{
    ofstream exportToFile;
    cout << "Writing File: " << fileName << endl;

    if (!append) {
        exportToFile.open(fileName);
        exportToFile << "picture ID, object ID, length, width, area, perimeter, depth, volume, posx, posY\n";
    } else {
        exportToFile.open(fileName, fstream::out | fstream::app);
    }
    
    // Loops through each piece and prints data to excel file
    for (int i = 0; i < pieces.size(); i++) {
        exportToFile << pieces[i].pictureID << ',' << pieces[i].objectID << ',' << pieces[i].length << ','
        << pieces[i].width << ',' << pieces[i].area << ',' << pieces[i].perimeter << ','
        << pieces[i].depth << ',' << pieces[i].volume << ',' << pieces[i].posX << ','
        << pieces[i].posY << endl;
    }
    
    exportToFile.close();
    cout << "File written" << endl;
}