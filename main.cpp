#include<gmp.h>
#include<iostream>

int main()
{
    mpz_t p, q, g, x, r, y, c1, c2, m, m1, m2;
    mpz_inits(p, q, g, x, y, r, c1, c2, m, m1, m2, NULL);

    mpz_set_str(p, "142282284370849912570755157394204324436151819483519736472506450277015504316442482302179809699487523684599695443823704505518649456607941635696764217139379994030733903992537187656778012130412386900614059558968181723449769746201162749723772776626562819525109211027308165765591038614679756134700779644992364964587", 10);
    mpz_set_str(q, "71141142185424956285377578697102162218075909741759868236253225138507752158221241151089904849743761842299847721911852252759324728303970817848382108569689997015366951996268593828389006065206193450307029779484090861724884873100581374861886388313281409762554605513654082882795519307339878067350389822496182482293", 10);
    mpz_set_str(g, "18975823582220682785993852932531386287495436262161294513590167289593267370356438398537208993434091045319278834737032652859606278917416498962717732710395691307744081603585840528342881011013680641549884113441794875506663788637433997099525153570862562610035270678141581938266029572392541284411362673338204862", 10);

    gmp_printf("p = %Zd\n", p);
    gmp_printf("q = %Zd\n", q);
    gmp_printf("g = %Zd\n", g);

    //n人とk人(n-1人)の判別。どこがnでどこがkか分ける

    int people, decode_people;
    people = 5;
    decode_people = people - 3;
    mpz_t key, is_slope[decode_people], is_x[people], is_f_x[people], tmp, pp, sum, sub;
    gmp_randstate_t stat;
    mpz_inits(tmp, pp, sub, NULL);
    gmp_randinit_mt(stat);

    mpz_init_set_str(key, "7", 10);
    mpz_init_set_str(sum, "0", 10);

    for (int i = 0; i < people; i++) // xに１から代入
    {
            mpz_init_set_str(is_x[i], "1", 10);
            mpz_mul_ui(is_x[i], is_x[i], i + 1);
    }
    for (int i = 0; i < decode_people; i++) // 傾きの乱数
    {
            mpz_init(is_slope[i]);
            mpz_urandomm(is_slope[i], stat, q);
    }
    for (int j = 0; j < people; j++) // f(x)の算出
    {
            mpz_init(is_f_x[j]);
            mpz_set(is_f_x[j], key);
            for (int i = 0; i < decode_people; i++)
            {
                    mpz_t is_tmp;
                    mpz_init_set_ui(is_tmp, 1);
                    mpz_pow_ui(is_tmp, is_x[j], i + 1);
                    mpz_mul(is_tmp, is_tmp, is_slope[i]); // i + 1
                    mpz_add(is_f_x[j], is_f_x[j], is_tmp);
                    mpz_clear(is_tmp);
            }
            mpz_mod(is_f_x[j], is_f_x[j], q);
    }

    // ラグランジュ補間

    for (int i = 0; i <= decode_people; i++)
    {
            mpz_set(pp, is_f_x[i]);
            for (int j = 0; j <= decode_people; j++)
            {
                    if (i != j)
                    {
                            mpz_neg(tmp, is_x[j]);
                            mpz_mul(pp, pp, tmp);
                            mpz_sub(sub, is_x[i], is_x[j]);
                            mpz_invert(tmp, sub, q);
                            mpz_mul(pp, pp, tmp);
                            mpz_mod(pp, pp, q);
                    }
            }
            mpz_add(sum, sum, pp);
            mpz_mod(sum, sum, q);
    }
    mpz_set(x, sum);
    mpz_powm(y, g, x, p);

    printf("\n");
    mpz_set_str(m, "12345", 10);
    gmp_printf("m1 = %Zd\n", m);

    mpz_urandomm(r, stat, q);
    mpz_powm(c1, g, r, p);
    mpz_powm(c2, y, r, p);
    mpz_mul(c2, c2, m);

    printf("\n");
    gmp_printf("c1 = %Zd\n", c1);
    gmp_printf("c2 = %Zd\n", c2);
    printf("\n");

    // 閾値復号

    mpz_t v[people];
    mpz_t(c1_s);
    mpz_init_set_str(c1_s, "1", 10);
    for(int i = 0; i < people; i++)
    {
            mpz_init_set_str(v[i], "1", 10);
            mpz_powm(v[i], c1, is_f_x[i], p);
    }
    for (int i = 0; i <= decode_people; i++)
    {
            mpz_set_str(pp, "1", 10);
            for (int j = 0; j <= decode_people; j++)
            {
                    if (i != j)
                    {
                            mpz_neg(tmp, is_x[j]); // -is_x[j]
                            mpz_mul(pp, pp, tmp);
                            mpz_sub(sub, is_x[i], is_x[j]);
                            mpz_invert(tmp, sub, q);
                            mpz_mul(pp, pp, tmp);
                            mpz_mod(pp, pp, q);
                    }
            }
            mpz_powm(tmp, v[i], pp, p);
            mpz_mul(c1_s, c1_s, tmp);
            mpz_mod(c1_s, c1_s, p);
    }
    mpz_invert(c1_s, c1_s, p);
    mpz_mul(m2, c2, c1_s);
    mpz_mod(m2, m2, p);
    gmp_printf("m2 = %Zd\n", m2);
    if(mpz_cmp(m, m2) == 0) {
        printf("success!!");
    } else {
        printf("failed......");
    }
}
