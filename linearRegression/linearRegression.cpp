#include "seal/seal.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>


#define debug_on
#ifdef debug_on
#define dbg(...) __va_args__
#else
#define dbg(...) 
#endif

using namespace seal;
using namespace std;

// 重新加密解密，并输出
void dbg_enc(Ciphertext& ct_1, double scale, Decryptor& decryptor, CKKSEncoder& encoder, Encryptor& encryptor) {
    
    Plaintext temp_p;
    decryptor.decrypt(ct_1, temp_p);
    vector<double> temp_v;
    encoder.decode(temp_p, temp_v);
    cout << "DBG" << endl;
    cout << temp_v.size() << endl;
    for (size_t i = 0; i < temp_v.size(); i++)
    {
            cout << temp_v[i] << " ";

    }

    encoder.encode(temp_v, scale, temp_p);
    encryptor.encrypt(temp_p, ct_1);
}
// 重新加密解密
void re_enc(Ciphertext& ct_1, double scale, Decryptor& decryptor, CKKSEncoder& encoder, Encryptor& encryptor) {
    
    Plaintext temp_p;
    decryptor.decrypt(ct_1, temp_p);
    vector<double> temp_v;
    encoder.decode(temp_p, temp_v);
    encoder.encode(temp_v, scale, temp_p);
    encryptor.encrypt(temp_p, ct_1);

}
//获取数据集的列数，注意，这里数据集规定为以逗号间隔的txt数据集，不要有标题栏
int getNumberOfColumns(const std::string& filename) {

    std::ifstream file(filename);
    std::string line;

    if (file.is_open()) {
        std::getline(file, line);
        file.close();

        std::istringstream iss(line);
        std::string token;
        int count = 0;

        while (std::getline(iss, token, ',')) {
            count++;
        }

        return count;
    }
    else {
        std::cout << "Failed to open file: " << filename << std::endl;
        return -1;
    }

}

int findNearestPowerOfTwo(int number) {
    
    if ((number & (number - 1)) == 0) {
        return number;
    }

    number--;
    for (int i = 1; i < sizeof(unsigned int) * 8; i <<= 1) {
        number |= (number >> i);
    }

    return number + 1;

}

int main() {

    // 打开文件
    std::ifstream file("train_set.txt");
    std::string line;

    // 定义feature和label容器
    std::vector<double> label;
    std::vector<double> feature;

    int col_num = getNumberOfColumns("train_set.txt") - 1; //特征的列数
    int col_num_of_vector = findNearestPowerOfTwo(col_num) * 4; //明文向量的列数，乘4是因为我的rotation的形式
    int row_num_of_vector = 4096 / col_num_of_vector; //明文向量的行数
    int row_num = 0; //数据集的行数

    while (std::getline(file, line)) {
        row_num ++;
        std::stringstream ss(line);
        std::string token;
        std::vector<double> row_data;

        // 提取每一行的数据
        while (std::getline(ss, token, ',')) {
            row_data.push_back(std::stod(token));
        }
        // 将第一列的label值push_back到label col_num + 1 次
        for (int i = 0; i < col_num + 1; i++) {
            label.push_back(row_data[0]);
        }
        // 为label每一行填充col_num_of_vector - col_num - 1 次0
        for (int i = 0; i < col_num_of_vector - col_num - 1; i++) {
            label.push_back((double)0);
        }
        // 将第二列到第col_num列的值push_back到feature
        for (int i = 1; i < col_num + 1; i++) {
            feature.push_back(row_data[i]);
        }

        //多push_back一次1，是为了加入偏置项
        feature.push_back(1);

        // 特征由 col_num + 1 扩展到 col_num_of_vector 列
        for (int i = 0; i < col_num_of_vector - col_num - 1; i++) {
            feature.push_back((double)0);
        }
    }

    // 扩展到4096 4096 - col_num_of_vector * row_num
    for (int i = 0; i < 4096 - col_num_of_vector * row_num; i++) {
        label.push_back((double)0);
        feature.push_back((double)0);
    }

    //// 输出检查
    //std::cout << "Label: ";
    //for (const auto& value : label) {
    //    std::cout << value << " ";
    //}
    //std::cout << std::endl;

    //std::cout << "Feature: ";
    //for (const auto& value : feature) {
    //    std::cout << value << " ";
    //}
    //std::cout << std::endl;


    //设置超参
    double lr = 0.6; // 学习率
    int epochs = 40; // 迭代次数
    double m = row_num; // 样本数量

    // 初始化权重向量、division(用以除样本数量)和learning rate
    //权重向量需要push_back col_num + 1次1，col_num_of_vector - col_num - 1次0，一共要row_num_of_vector次
    vector<double> weight;
    for (int i = 0; i < row_num_of_vector; i++) {
        
        for (int i = 0; i < col_num + 1; i++) {
            weight.push_back(1);
        }
        for (int i = 0; i < col_num_of_vector - col_num - 1; i++) {
            weight.push_back(0);
        }
    }
    vector<double> division(4096, (double)1 / m);
    vector<double> lr_v(4096, lr);

    // 前置设置
    EncryptionParameters parms(scheme_type::ckks);
    size_t poly_modulus_degree = 8192;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 30, 30, 30, 60 }));
    double scale = pow(2.0, 30);
    SEALContext context(parms);
    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);
    //keygen.create_galois_keys(setps，gal_keys) /////////////////////////////////////////////////////////注意这里，要用这种方式
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    CKKSEncoder encoder(context);
    //size_t slot_count = encoder.slot_count();
    //cout << "Number of slots: " << slot_count << endl;


    // 加密feature
    Plaintext x_plain1;
    encoder.encode(feature, scale, x_plain1);
    Ciphertext ct_x;
    encryptor.encrypt(x_plain1, ct_x);

    // 加密label
    Plaintext x_plain2;
    encoder.encode(label, scale, x_plain2);
    Ciphertext ct_y;
    encryptor.encrypt(x_plain2, ct_y);


    // 加密weight
    Plaintext x_plain3;
    encoder.encode(weight, scale, x_plain3);
    Ciphertext ct_b;
    encryptor.encrypt(x_plain3, ct_b);

    // 加密数据组数
    Plaintext x_plain4;
    encoder.encode(division, scale, x_plain4);
    Ciphertext ct_division;
    encryptor.encrypt(x_plain4, ct_division);

    // 加密学习率
    Plaintext x_plain5;
    encoder.encode(lr_v, scale, x_plain5);
    Ciphertext ct_lr;
    encryptor.encrypt(x_plain5, ct_lr);

    /*dbg_enc(ct_b, scale, decryptor, encoder, encryptor)*/;

    for (int j = 1; j <= epochs; j++) {
        // 计算h0i = = b * x
        Ciphertext ct_h0i;
        evaluator.multiply(ct_b, ct_x, ct_h0i);
        evaluator.relinearize_inplace(ct_h0i, relin_keys);
        evaluator.rescale_to_next_inplace(ct_h0i);
        //dbg_enc(ct_h0i, scale, decryptor, encoder, encryptor);

        //计算sum h0i
        Ciphertext ct_h0 = ct_h0i;
        Ciphertext temp1;
        Ciphertext temp2;
        for (int i = 1; i < col_num + 1; i++)
        {
            evaluator.rotate_vector(ct_h0i, -i, gal_keys, temp1);
            evaluator.rotate_vector(ct_h0i, i, gal_keys, temp2);
            evaluator.add_inplace(ct_h0, temp1);
            evaluator.add_inplace(ct_h0, temp2);
        }
        re_enc(ct_h0, scale, decryptor, encoder, encryptor);
        //dbg_enc(ct_h0, scale, decryptor, encoder, encryptor);

        //计算y - h0
        Ciphertext ct_y_sub_h0;
        evaluator.sub(ct_y, ct_h0, ct_y_sub_h0);
        //dbg_enc(ct_y_sub_h0, scale, decryptor, encoder, encryptor);

        //计算x * (y - h0)
        Ciphertext ct_x_multi_y_sub_h0;
        evaluator.multiply(ct_x, ct_y_sub_h0, ct_x_multi_y_sub_h0);
        evaluator.relinearize_inplace(ct_x_multi_y_sub_h0, relin_keys);
        evaluator.rescale_to_next_inplace(ct_x_multi_y_sub_h0);
        //dbg_enc(ct_x_multi_y_sub_h0, scale, decryptor, encoder, encryptor);

        //计算sum上边那个
        Ciphertext ct_sum = ct_x_multi_y_sub_h0;
        Ciphertext temp3;
        Ciphertext temp4;
        int rotate = col_num_of_vector;
        for (int i = 1; i < m; i++)
        {
            evaluator.rotate_vector(ct_x_multi_y_sub_h0, rotate * i, gal_keys, temp3);
            evaluator.rotate_vector(ct_x_multi_y_sub_h0, - rotate * i, gal_keys, temp4);
            evaluator.add_inplace(ct_sum, temp3);
            evaluator.add_inplace(ct_sum, temp4);
        }
        re_enc(ct_sum, scale, decryptor, encoder, encryptor);
        //dbg_enc(ct_sum, scale, decryptor, encoder, encryptor);

        //计算sum/m
        Ciphertext ct_sum_div_m;
        evaluator.multiply(ct_sum, ct_division, ct_sum_div_m);
        evaluator.relinearize_inplace(ct_sum_div_m, relin_keys);
        evaluator.rescale_to_next_inplace(ct_sum_div_m);
        re_enc(ct_sum_div_m, scale, decryptor, encoder, encryptor);
        //dbg_enc(ct_sum_div_m, scale, decryptor, encoder, encryptor);

        //上面那个乘lr
        Ciphertext ct_multi_lr;
        evaluator.multiply(ct_sum_div_m, ct_lr, ct_multi_lr);
        evaluator.relinearize_inplace(ct_multi_lr, relin_keys);
        evaluator.rescale_to_next_inplace(ct_multi_lr);
        re_enc(ct_multi_lr, scale, decryptor, encoder, encryptor);
        //dbg_enc(ct_multi_lr, scale, decryptor, encoder, encryptor);

        //改变权重
        evaluator.add_inplace(ct_b, ct_multi_lr);

        std::cout << "第" << j << "次迭代完成" << endl;
    }

    //得到的权重（前四个feature的weight加第五个bias）
    Plaintext temp_p;
    decryptor.decrypt(ct_b, temp_p);
    vector<double> temp_v;
    encoder.decode(temp_p, temp_v);

    cout << "权重:" << temp_v[0] <<" " << temp_v[1] << " " << temp_v[2] << "" << temp_v[3] << " " << temp_v[4] << endl;

    //利用测试集进行预测
    std::string filename = "test_set.txt";
    std::ifstream file1(filename);
    std::string line1;

    double loss = 0;
    int test_num = 0;
    if (file1.is_open()) {
        while (std::getline(file1, line1)) {
            test_num++;
            std::vector<double> values;
            std::istringstream iss(line1);
            std::string token;

            while (std::getline(iss, token, ',')) {
                double value = std::stod(token);
                values.push_back(value);
            }

            if (values.size() == 5) {
                values[1] *= temp_v[0];
                values[2] *= temp_v[1];
                values[3] *= temp_v[2];
                values[4] *= temp_v[3];

                double result = values[1] + values[2] + values[3] + values[4] + temp_v[4];
                std::cout << std::fixed << std::setprecision(4);
                std::cout << "实际值" << values[0] << "预测值" << result << endl;
                loss += sqrt((values[0] - result) * (values[0] - result));
            }
            
        }

        std::cout << "平均误差" << loss / test_num;
        file1.close();
    }
}
