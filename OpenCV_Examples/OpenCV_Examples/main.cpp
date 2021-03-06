/* This is sample from the OpenCV book. The copyright notice is below */

/* *************** License:**************************
 Oct. 3, 2008
 Right to use this code in any way you want without warranty, support or any guarantee of it working.
 
 BOOK: It would be nice if you cited it:
 Learning OpenCV: Computer Vision with the OpenCV Library
 by Gary Bradski and Adrian Kaehler
 Published by O'Reilly Media, October 3, 2008
 
 AVAILABLE AT:
 http://www.amazon.com/Learning-OpenCV-Computer-Vision-Library/dp/0596516134
 Or: http://oreilly.com/catalog/9780596516130/
 ISBN-10: 0596516134 or: ISBN-13: 978-0596516130
 
 OPENCV WEBSITES:
 Homepage:      http://opencv.org
 Online docs:   http://docs.opencv.org
 Q&A forum:     http://answers.opencv.org
 Issue tracker: http://code.opencv.org
 GitHub:        https://github.com/Itseez/opencv/
 ************************************************** */
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/contrib/contrib.hpp"

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

using namespace cv;
using namespace std;

static void saveXYZ(const char* filename, const Mat& mat)
{
    const double max_z = 1.0e4;
    FILE* fp = fopen(filename, "wt");
    for(int y = 0; y < mat.rows; y++)
    {
        for(int x = 0; x < mat.cols; x++)
        {
            Vec3f point = mat.at<Vec3f>(y, x);
            if(fabs(point[2] - max_z) < FLT_EPSILON || fabs(point[2]) > max_z) continue;
            fprintf(fp, "%f %f %f\n", point[0], point[1], point[2]);
        }
    }
    fclose(fp);
}

static bool readStringList( const string& filename, vector<string>& l )
{
    l.resize(0);
    FileStorage fs(filename, FileStorage::READ);
    if( !fs.isOpened() )
        return false;
    
    FileNode n = fs.getFirstTopLevelNode();
    if( n.type() != FileNode::SEQ )
        return false;
    
    FileNodeIterator it = n.begin(), it_end = n.end();
    for( ; it != it_end; ++it )
        l.push_back((string)*it);
    
    return true;
}

static void stereoCalCap(Size boardsize)
{
    VideoCapture capR(0), capL(1);
    Mat right, left, grayR, grayL;
    Mat cornersR, cornersL;
    
    int k, count = 0, numpairs = 6;
    bool foundR = false, foundL = false;
    
    string fileR, fileL;
    string path = "/Users/vegas_bballer/Documents/Senior_Project/Testing/images/image";
    string ext = ".jpg";
    string r = "Right";
    string l = "Left";
    
    capR.set(CV_CAP_PROP_FRAME_HEIGHT, 720);
    capR.set(CV_CAP_PROP_FRAME_WIDTH, 1280);
    
    capL.set(CV_CAP_PROP_FRAME_HEIGHT, 720);
    capL.set(CV_CAP_PROP_FRAME_WIDTH, 1280);
    
    while(count < numpairs)
    {
        //capR.read(right);
        //capL.read(left);
        
        capR >> right;
        capL >> left;
        
        //resize(right, right, Size(800, 700));
        //resize(left, left, Size(800, 700));
        
        cvtColor(right, grayR, CV_RGB2GRAY);
        cvtColor(left, grayL, CV_RGB2GRAY);
        
        foundR = findChessboardCorners(right, boardsize, cornersR, CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_NORMALIZE_IMAGE);
        foundL = findChessboardCorners(left, boardsize, cornersL, CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_NORMALIZE_IMAGE);
        
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

static void stereoCalib(const vector<string>& imagelist, Size boardSize, bool useCalibrated=true, bool showRectified=true)
{
    if( imagelist.size() % 2 != 0 )
    {
        cout << "Error: the image list contains odd (non-even) number of elements\n";
        return;
    }
    
    bool displayCorners = false;//true;
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
                if( found )
                {
                    if( scale > 1 )
                    {
                        Mat cornersMat(corners);
                        cornersMat *= 1./scale;
                    }
                    break;
                }
            }
            if( displayCorners )
            {
                cout << filename << endl;
                Mat cimg, cimg1;
                cvtColor(img, cimg, COLOR_GRAY2BGR);
                drawChessboardCorners(cimg, boardSize, corners, found);
                double sf = 640./MAX(img.rows, img.cols);
                resize(cimg, cimg1, Size(), sf, sf);
                imshow("corners", cimg1);
                //char c = (char)waitKey(500);
                //if( c == 27 || c == 'q' || c == 'Q' ) //Allow ESC to quit
                //  exit(-1);
                imwrite("/Users/vegas_bballer/Documents/Senior_Project/Testing/images/corners.jpg", cimg1);
                waitKey(0);
            }
            else
                putchar('.');
            if( !found )
                break;
            cornerSubPix(img, corners, Size(11,11), Size(-1,-1),
                         TermCriteria(CV_TERMCRIT_ITER+CV_TERMCRIT_EPS,
                                      30, 0.01));
        }
        if( k == 2 )
        {
            goodImageList.push_back(imagelist[i*2]);
            goodImageList.push_back(imagelist[i*2+1]);
            j++;
        }
    }
    cout << j << " pairs have been successfully detected.\n";
    nimages = j;
    if( nimages < 2 )
    {
        cout << "Error: too little pairs to run the calibration\n";
        return;
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
    FileStorage fs("/Users/vegas_bballer/Documents/Senior_Project/Testing/intrinsics.yml", CV_STORAGE_WRITE);
    if( fs.isOpened() )
    {
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
    
    fs.open("/Users/vegas_bballer/Documents/Senior_Project/Testing/extrinsics.yml", CV_STORAGE_WRITE);
    if( fs.isOpened() )
    {
        fs << "R" << R << "T" << T << "R1" << R1 << "R2" << R2 << "P1" << P1 << "P2" << P2 << "Q" << Q;
        fs.release();
    }
    else
        cout << "Error: can not save the extrinsic parameters\n";
    
    // OpenCV can handle left-right
    // or up-down camera arrangements
    bool isVerticalStereo = fabs(P2.at<double>(1, 3)) > fabs(P2.at<double>(0, 3));
    
    // COMPUTE AND DISPLAY RECTIFICATION
    if( !showRectified )
        return;
    
    Mat rmap[2][2];
    // IF BY CALIBRATED (BOUGUET'S METHOD)
    if( useCalibrated )
    {
        // we already computed everything
        cout << "Bouguet's Method\n";
    }
    // OR ELSE HARTLEY'S METHOD
    else
        // use intrinsic parameters of each camera, but
        // compute the rectification transformation directly
        // from the fundamental matrix
    {
        vector<Point2f> allimgpt[2];
        for( k = 0; k < 2; k++ )
        {
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
    }
    
    //Precompute maps for cv::remap()
    initUndistortRectifyMap(cameraMatrix[0], distCoeffs[0], R1, P1, imageSize, CV_16SC2, rmap[0][0], rmap[0][1]);
    initUndistortRectifyMap(cameraMatrix[1], distCoeffs[1], R2, P2, imageSize, CV_16SC2, rmap[1][0], rmap[1][1]);
    
    Mat canvas;
    double sf;
    int w, h;
    if( !isVerticalStereo )
    {
        sf = 600./MAX(imageSize.width, imageSize.height);
        w = cvRound(imageSize.width*sf);
        h = cvRound(imageSize.height*sf);
        canvas.create(h, w*2, CV_8UC3);
    }
    else
    {
        sf = 300./MAX(imageSize.width, imageSize.height);
        w = cvRound(imageSize.width*sf);
        h = cvRound(imageSize.height*sf);
        canvas.create(h*2, w, CV_8UC3);
    }
    
    for( i = 0; i < nimages; i++ )
    {
        for( k = 0; k < 2; k++ ) // left then right
        {
            Mat img = imread(goodImageList[i*2+k], 0), rimg, cimg;
            
            imshow("IMG", img);
            waitKey(0);
            
            remap(img, rimg, rmap[k][0], rmap[k][1], CV_INTER_LINEAR);
            cvtColor(rimg, cimg, COLOR_GRAY2BGR);
            Mat canvasPart = canvas(Rect(w*k, 0, w, h));
            resize(cimg, canvasPart, canvasPart.size(), 0, 0, CV_INTER_AREA);
            //if( useCalibrated )
            {
                
                //Rect vroi(cvRound(validRoi[k].x*sf), cvRound(validRoi[k].y*sf),
                  //        cvRound(validRoi[k].width*sf), cvRound(validRoi[k].height*sf));
                //rectangle(canvasPart, vroi, Scalar(0,0,255), 3, 8);
            }
            cout << i << " " << k << endl;
        }
        
        //if( !isVerticalStereo )
            for( j = 0; j < canvas.rows; j += 16 )
                line(canvas, Point(0, j), Point(canvas.cols, j), Scalar(0, 255, 0), 1, 8);
        //else
            //for( j = 0; j < canvas.cols; j += 16 )
              //  line(canvas, Point(j, 0), Point(j, canvas.rows), Scalar(0, 255, 0), 1, 8);
        
        imshow("Rectified", canvas);
        
        
        
        //string file = "/Users/vegas_bballer/Documents/Senior_Project/Testing/rectified" + to_string(i + 1) + ".jpg";
        //imwrite(file, canvas);
        
        
        char c = (char)waitKey();
        if( c == 27 || c == 'q' || c == 'Q' )
            break;
    }
    destroyAllWindows();
}

static void stereoCap()
{
    VideoCapture capR(0), capL(1);
    Mat right, left;
    
    string fileR, fileL;
    string path = "/Users/vegas_bballer/Documents/Senior_Project/Testing/images/capture";
    string ext = ".jpg";
    string r = "Right";
    string l = "Left";
    int count = 0, numpics = 2, k;
    
    capR.set(CV_CAP_PROP_FRAME_HEIGHT, 720);
    capR.set(CV_CAP_PROP_FRAME_WIDTH, 1280);
    
    capL.set(CV_CAP_PROP_FRAME_HEIGHT, 720);
    capL.set(CV_CAP_PROP_FRAME_WIDTH, 1280);
    
    while(count < numpics)
    {
        //capR.read(right);
        //capL.read(left);
        
        capR >> right;
        capL >> left;
        
        //resize(right, right, Size(800, 700));
        //resize(left, left, Size(800, 700));
        
        imshow("Capture Right", right);
        imshow("Capture Left", left);
        
        k = waitKey(10);
        
        if (k == 13 ) // return key
        {
            cout << "Saving images.\n";
            fileR = path + r + to_string(count + 1) + ext;
            fileL = path + l + to_string(count + 1) + ext;
            cout << fileL << endl << fileR << endl;
            
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

int stereoMatch()
{
    const char* imgL_filename = "/Users/vegas_bballer/Documents/Senior_Project/Testing/images/captureLeft1.jpg";
    const char* imgR_filename = "/Users/vegas_bballer/Documents/Senior_Project/Testing/images/captureRight1.jpg";
    const char* intrinsic_filename = "/Users/vegas_bballer/Documents/Senior_Project/Testing/intrinsics.yml";
    const char* extrinsic_filename = "/Users/vegas_bballer/Documents/Senior_Project/Testing/extrinsics.yml";
    const char* disparity_filename = "/Users/vegas_bballer/Documents/Senior_Project/Testing/disparity.jpg";
    const char* point_cloud_filename = "/Users/vegas_bballer/Documents/Senior_Project/Testing/pointcloud.jpg";
    
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
    Mat imgL = imread(imgL_filename, color_mode);
    Mat imgR = imread(imgR_filename, color_mode);
    
    Size img_size = imgL.size();
    Rect roi1, roi2;
    Mat Q;
    
    if( intrinsic_filename )
    {
        // reading intrinsic parameters
        FileStorage fs(intrinsic_filename, CV_STORAGE_READ);
        if(!fs.isOpened())
        {
            printf("Failed to open file %s\n", intrinsic_filename);
            return -1;
        }
        
        Mat M1, D1, M2, D2;
        fs["M1"] >> M1;
        fs["D1"] >> D1;
        fs["M2"] >> M2;
        fs["D2"] >> D2;
        
        M1 *= scale;
        M2 *= scale;
        
        fs.open(extrinsic_filename, CV_STORAGE_READ);
        if(!fs.isOpened())
        {
            printf("Failed to open file %s\n", extrinsic_filename);
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
    
    //numberOfDisparities = numberOfDisparities > 0 ? numberOfDisparities : ((img_size.width/8) + 15) & -16;
    
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
    
    //sgbm.preFilterCap = 63;
    //sgbm.SADWindowSize = SADWindowSize > 0 ? SADWindowSize : 3;
    
    //int cn = imgL.channels();
    
    //sgbm.P1 = 8*cn*sgbm.SADWindowSize*sgbm.SADWindowSize;
    //sgbm.P2 = 32*cn*sgbm.SADWindowSize*sgbm.SADWindowSize;
    //sgbm.minDisparity = 0;
    //sgbm.numberOfDisparities = numberOfDisparities;
    //sgbm.uniquenessRatio = 30;
    //sgbm.speckleWindowSize = bm.state->speckleWindowSize;
    //sgbm.speckleRange = bm.state->speckleRange;
    //sgbm.disp12MaxDiff = 1;
    //sgbm.fullDP = alg == STEREO_SGBM;
    
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
    
    sgbm.SADWindowSize = 5; //3
    sgbm.numberOfDisparities = numberOfDisparities; //11*16
    sgbm.preFilterCap = 80; //15
    sgbm.minDisparity = 0; //5
    sgbm.uniquenessRatio = 8; //20
    sgbm.speckleWindowSize = 1000; //30
    sgbm.speckleRange = 4; //2
    sgbm.disp12MaxDiff = 0; //100
    sgbm.fullDP = false; //true
    sgbm.P1 = 3000; //200
    sgbm.P2 = 10000; //3000
    
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
    
    //int64 t = getTickCount();
    if( alg == STEREO_BM )
        bm(imgL, imgR, disp);
    else if( alg == STEREO_VAR ) {
        var(imgL, imgR, disp);
    }
    else if( alg == STEREO_SGBM || alg == STEREO_HH )
        sgbm(imgL, imgR, disp);
    
    //t = getTickCount() - t;
    //printf("Time elapsed: %fms\n", t*1000/getTickFrequency());
    
    //disp = dispp.colRange(numberOfDisparities, img1p.cols);
    if( alg != STEREO_VAR )
    {
        disp.convertTo(disp8, CV_8U, 255/(numberOfDisparities*16.));
        cout << 255/(numberOfDisparities * 16.) << endl;
    }
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
        
        printf("press any key to continue...");
        fflush(stdout);
        waitKey();
        printf("\n");
    }
    
    if(disparity_filename)
        imwrite(disparity_filename, disp8);
    
    if(point_cloud_filename)
    {
        printf("storing the point cloud...");
        fflush(stdout);
        Mat xyz;
        reprojectImageTo3D(disp, xyz, Q, true);
        saveXYZ(point_cloud_filename, xyz);
        printf("\n");
    }
    destroyAllWindows();
    return 0;
}


int main()
{
    Size boardSize = Size(9, 6);
    string imagelistfn = "/Users/vegas_bballer/Documents/Senior_Project/Testing/stereo_calib.xml";
    bool showRectified = true;
    
    vector<string> imagelist;
    bool ok = readStringList(imagelistfn, imagelist);
    
    if (!ok || imagelist.empty())
    {
        cout << "can not open " << imagelistfn << " or the string list is empty" << endl;
        return -1;
    }
    
    //stereoCalCap(boardSize);
    //stereoCalib(imagelist, boardSize, false, showRectified);
    //stereoCap();    // Get pictures to view disparity map
    stereoMatch();  // Calculate a disparty map
    //while(1)
    {
      //  return 0;
    }
    bool done = false;
    int k = 0;
    StereoSGBM sgbm;
    Mat disp, disp8;
    Mat imgL = imread("/Users/vegas_bballer/Documents/Senior_Project/Testing/images/captureLeft1.jpg");
    Mat imgR = imread("/Users/vegas_bballer/Documents/Senior_Project/Testing/images/captureRight1.jpg");
    string intrinsicFile = "/Users/vegas_bballer/Documents/Senior_Project/Testing/intrinsics.yml";
    string extrinsicFile = "/Users/vegas_bballer/Documents/Senior_Project/Testing/extrinsics.yml";
    //imshow("Right Image", imgR);
    //imshow("Left Image", imgL);
    
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
        
        //imshow("IMG LR",imgLr);
        //waitKey(0);
        //imshow("IMG L",imgL);
        //waitKey(0);
        //imshow("IMG RR",imgRr);
        //waitKey(0);
        //imshow("IMG R",imgR);
        //waitKey(0);
        
        imgL = imgLr;
        imgR = imgRr;
    }
    else
        cout << "File not opened\n";
    
    int scalar = 16;
    const double SADWindowSize_max = 15;
    const double numberOfDisparities_max = 20;
    const double preFilterCap_max = 100;
    const double uniquenessRatio_max = 50;
    const double speckleWindowSize_max = 200;
    const double speckleRange_max = 5;
    const double disp12MaxDiff_max = 300;
    const double P1_max = 4000;
    const double P2_max = 12000;
    
    namedWindow("Stereo Tune");
    
    sgbm.SADWindowSize = 5; //3
    sgbm.numberOfDisparities = scalar * 16; //11*16
    sgbm.preFilterCap = 6; //15
    sgbm.minDisparity = 0; //5
    sgbm.uniquenessRatio = .08; //20
    sgbm.speckleWindowSize = 30; //30
    sgbm.speckleRange = 4; //2
    sgbm.disp12MaxDiff = 0; //100
    sgbm.fullDP = false; //true
    sgbm.P1 = 3000; //200
    sgbm.P2 = 12000; //3000
    
    createTrackbar("SAD Window Size",       "Stereo Tune", &sgbm.SADWindowSize,         SADWindowSize_max);
    //createTrackbar("Number of Disparities", "Stereo Tune", &scalar,                     numberOfDisparities_max);
    //sgbm.numberOfDisparities = scalar * 16;
    createTrackbar("Pre Filter Cap",        "Stereo Tune", &sgbm.preFilterCap,          preFilterCap_max);
    createTrackbar("Uniqueness Ratio",      "Stereo Tune", &sgbm.uniquenessRatio,       uniquenessRatio_max);
    createTrackbar("Speckle Window Size",   "Stereo Tune", &sgbm.speckleWindowSize,     speckleWindowSize_max);
    createTrackbar("Speckle Range",         "Stereo Tune", &sgbm.speckleRange,          speckleRange_max);
    createTrackbar("disp12MaxDiff_Max",     "Stereo Tune", &sgbm.disp12MaxDiff,         disp12MaxDiff_max);
    createTrackbar("P1",                    "Stereo Tune", &sgbm.P1,                    P1_max);
    createTrackbar("P2",                    "Stereo Tune", &sgbm.P2,                    P2_max);
    
    //sgbm.preFilterCap = 63;
    //sgbm.SADWindowSize = SADWindowSize > 0 ? SADWindowSize : 3; //19
    //int cn = imgL.channels();
    //sgbm.P1 = 8*cn*sgbm.SADWindowSize*sgbm.SADWindowSize; // 2888
    //sgbm.P2 = 32*cn*sgbm.SADWindowSize*sgbm.SADWindowSize; // 11552
    //sgbm.minDisparity = 0;
    //sgbm.numberOfDisparities = numberOfDisparities; // 176
    //sgbm.uniquenessRatio = 30; //or 10
    //sgbm.speckleWindowSize = bm.state->speckleWindowSize; // 0
    //sgbm.speckleRange = bm.state->speckleRange; // 1
    //sgbm.disp12MaxDiff = 1;
    //sgbm.fullDP = alg == STEREO_HH;
    
    while (!done)
    {
        cout << "while loop\n";
        sgbm(imgL, imgR, disp);
        imshow("Disp",disp);
        disp.convertTo(disp8, CV_8U, 0.2);
        
        imshow("Left", imgL);
        imshow("Right", imgR);
        namedWindow("Disparity");
        imshow("Disparity", disp8);
        imwrite("/Users/vegas_bballer/Documents/Senior_Project/Testing/Output/disparity.jpg", disp8);
        imwrite("/Users/vegas_bballer/Documents/Senior_Project/Testing/Output/rightImage.jpg", imgR);
        imwrite("/Users/vegas_bballer/Documents/Senior_Project/Testing/Output/leftImage.jpg", imgL);
        
        waitKey(0);
        //k = waitKey(10);
        //if (k == 13)
            //done = true;
        
    }
    return 0;
}








