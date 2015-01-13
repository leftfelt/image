#include <gtest/gtest.h>
#include <Image.hpp>
#include <ImageIO.hpp>

TEST(ImageTest, PutAndGet){
    Image image1(2,2,std::vector<unsigned char>{2,4,4,2});
    Image image2(2,2,std::vector<unsigned char>{1,2,2,1});

    ASSERT_EQ(Pixel(2), image1.Get(0,0));
    ASSERT_EQ(Pixel(4), image1.Get(1,0));
    ASSERT_EQ(Pixel(4), image1.Get(0,1));
    ASSERT_EQ(Pixel(2), image1.Get(1,1));

    image1.Put(0,1,Pixel(255));

    ASSERT_EQ(Pixel(255), image1.Get(0,1));
}

TEST(ImageTest, OperatorImageAndImage){
    Image image1(2,2,std::vector<unsigned char>{2,4,4,2});
    Image image2(2,2,std::vector<unsigned char>{1,2,2,1});

    ASSERT_EQ(Image(2,2,std::vector<unsigned char>{3,6,6,3}), image1 + image2);
    ASSERT_EQ(Image(2,2,std::vector<unsigned char>{1,2,2,1}), image1 - image2);
    ASSERT_EQ(Image(2,2,std::vector<unsigned char>{2,8,8,2}), image1 * image2);
    ASSERT_EQ(Image(2,2,std::vector<unsigned char>{2,2,2,2}), image1 / image2);
}

TEST(ImageTest, OperatorImageAndPixel){
    Image image1(2,2,std::vector<unsigned char>{2,4,4,2});
    Image image2(2,2,std::vector<unsigned char>{1,2,2,1});

    ASSERT_EQ(Image(2,2,std::vector<unsigned char>{3,5,5,3}), image1 + Pixel(1));
    ASSERT_EQ(Image(2,2,std::vector<unsigned char>{1,3,3,1}), image1 - Pixel(1));
    ASSERT_EQ(Image(2,2,std::vector<unsigned char>{4,8,8,4}), image1 * Pixel(2));
    ASSERT_EQ(Image(2,2,std::vector<unsigned char>{1,2,2,1}), image1 / Pixel(2));
}


TEST(ImageTest, PasteTest){
    Image image1(2,2,std::vector<unsigned char>{2,4,4,2});
    Image image2(2,2,std::vector<unsigned char>{1,2,2,1});

    image1.Paste(0, 1, image2);
    ASSERT_EQ(Image(2,2,std::vector<unsigned char>{2,1,1,2}), image1); //image2の下にはみ出した部分は折り返してimage1の上部にPasteされる
}


TEST(ImageTest, cut){
    Image image1(2,2,std::vector<unsigned char>{2,4,4,2});

    ASSERT_EQ(Image(2,1,std::vector<unsigned char>{4,2}), image1.Cut(0,1,2,1));
}

TEST(ImageTest, Rotate){
    Image image1(2,2,std::vector<unsigned char>{0,255,255,0});
    Image image2(3,3,std::vector<unsigned char>{
            1,1,1,
            2,3,2,
            1,1,1
    });
    Image image3(4,4,std::vector<unsigned char>{
            1,1,1,1,
            2,3,3,2,
            2,3,3,2,
            1,1,1,1
    });

    ASSERT_EQ(Image(2,2,std::vector<unsigned char>{0,0,255,0}), image1.Rotate(90));
    ASSERT_EQ(Image(3,3,std::vector<unsigned char>{
            1,1,1,
            2,3,2,
            1,1,1
    }), image2.Rotate(180));
    ASSERT_EQ(Image(4,4,std::vector<unsigned char>{
            0,0,0,0,
            0,1,1,1,
            0,2,3,3,
            0,2,3,3,
    }), image3.Rotate(180));
    ASSERT_EQ(Image(4,4,std::vector<unsigned char>{
            0,0,0,0,
            1,2,2,1,
            1,3,3,1,
            1,3,3,1
    }), image3.Rotate(90));
}


TEST(ImageTest, Clear){
    Image image1(4,4,std::vector<unsigned char>{
            1,1,1,1,
            2,3,3,2,
            2,3,3,2,
            1,1,1,1
    });
    Image image2(2,2,std::vector<unsigned char>{
            0,255,
            255,0
    });
    image2.Clear(Pixel(0));
    image1.Paste(1,1,image2);
    ASSERT_EQ(Image(4,4,std::vector<unsigned char>{
            1,1,1,1,
            2,3,255,2,
            2,255,3,2,
            1,1,1,1
    }), image1);
}

TEST(ImageTest, Size){
    Image image1(4,4,std::vector<unsigned char>{
            1,1,1,1,
            2,3,3,2,
            2,3,3,2,
            1,1,1,1
    });
    image1.Size(3,3);
    ASSERT_EQ(Image(3,3,std::vector<unsigned char>{
            1,1,1,
            2,3,2,
            1,1,1
    }), image1);
    image1.Size(4,4);
    ASSERT_EQ(Image(4,4,std::vector<unsigned char>{
            1,1,1,1,
            2,3,2,2,
            1,1,1,1,
            1,1,1,1
    }), image1);
    image1.Size(2,2);
    ASSERT_EQ(Image(2,2,std::vector<unsigned char>{
            1,1,
            1,1
    }), image1);
}