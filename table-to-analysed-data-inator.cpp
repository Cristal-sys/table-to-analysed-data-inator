#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<cmath>
#include<vector>

struct Variable {

std::string name="0";//name

long double value=0;

long double sigma=0;

long double standard_deviation=0;

long double standard_deviation_mean=0;

};

long double expression_evaluator(std::string);
long double operacao(long double, char, long double);
std::string return_operand_string(std::stringstream&);
long double function_to_double(std::string f, std::vector<Variable> parameters);
long double partial(std::string f, std::string variable, std::vector<Variable> parameters);
void read_parameters(std::fstream& file, std::vector<Variable>&parameters);
long double uncertainty(std::string f, std::vector<Variable> parameters);
void read_variable(std::fstream& file, std::vector<Variable>& parameters);
void read_txt_table(std::fstream& file, std::vector<Variable>& parameters);
void read_function(std::fstream& file, std::vector<std::string>& functions);
void click_anything();
void receptionist(std::string& filename, std::string& command);
std::string get_edited_table_name(std::string filename);
void write(std::fstream& output, std::string command, std::vector<Variable> parameters, std::vector<std::string> functions);
long double angular_right(std::vector<Variable> parameters, long double n);
long double linear_right(std::vector<Variable> parameters, long double n);
long double angular_sigma_right(std::vector<Variable> parameters, long double n);
long double linear_sigma_right(std::vector<Variable> parameters, long double n);
long double error(std::vector<Variable> parameters);
long double a_2(std::vector<Variable> parameters, long double n);
long double a_1(std::vector<Variable> parameters, long double n);
long double a_0(std::vector<Variable> parameters, long double n);
long double a_2_sigma(std::vector<Variable> parameters, long double n);
long double a_1_sigma(std::vector<Variable> parameters, long double n);
long double a_0_sigma(std::vector<Variable> parameters, long double n);
long double error_parabola(std::vector<Variable> parameters);
bool is_operator_or_parenthesis(char);
void read_operand(std::string& function, std::string& c, int& i);
void remove_commentary(std::fstream& filename, std::string aux);
void restore_commentary(std::fstream& filename, std::string& backup);

const std::string euler="2.718281828459045235360287471352";
const std::string pi="3.141592653589793238462643383279";

int main(){
std::string filename, command;
std::vector<Variable> parameters;
std::vector<std::string> functions;

receptionist(filename, command);
std::fstream file(filename, std::ios::in);
while(!file.is_open()){
    std::cout<<"\nfile does not exist. ";
    receptionist(filename, command);
    file.close();
    file.open(filename, std::ios::in);
}

std::string final_file_name=get_edited_table_name(filename);
std::fstream output(final_file_name, std::ios::out);

read_txt_table(file, parameters);
read_function(file, functions);

write(output, command, parameters, functions);

file.close();
output.close();
system(final_file_name.c_str());
click_anything();

return 0;
}

long double expression_evaluator(std::string expression)
{
    std::stringstream expressionstream{expression};
    long double operand1, operand2;
    char operador;
        if(expressionstream.peek()=='(')
        {
            operand1=expression_evaluator(return_operand_string(expressionstream));
        } else
        {
            expressionstream>>operand1;
        }

    expressionstream>>operador;

        if(expressionstream.peek()=='(')
        {
            operand2=expression_evaluator(return_operand_string(expressionstream));
        } else
        {
            expressionstream>>operand2;
        }

  return operacao(operand1, operador, operand2);
}

long double operacao(long double operand1, char operador, long double operand2)
{
    if(operador=='+')
    {
        return operand1+operand2;
    }

    if(operador=='-')
    {
        return operand1-operand2;
    }

    if(operador=='*')
    {
        return operand1*operand2;
    }

    if(operador=='/')
    {
        return operand1/operand2;
    }

    if(operador=='^')
    {
        return std::pow(operand1, operand2);
    }

    if(operador=='S')
    {
	if(operand1==0){
		return std::asin(operand2);
	}
        return std::sin(operand2);
    }

    if(operador=='C')
    {
	if(operand1==0){
		return std::acos(operand2);
	}
        return std::cos(operand2);
    }
    if(operador=='T')
    {
	if(operand1==0){
		return std::atan(operand2);
	}
        return std::tan(operand2);
    }

    if(operador=='L')
    {
        return std::log(operand2);
    }

    return 0;
}

std::string return_operand_string(std::stringstream& expressionstream)
{
    std::string operand;
            char c;
            expressionstream>>c;
            int x=1;
            while(true)
            {
                if(x==1 && expressionstream.peek()==')')
                {
                    expressionstream>>c;
                    break;
                }
                expressionstream>>c;
                if(c=='(')
                {
                    ++x;
                }
                if(c==')')
                {
                    --x;
                }
                operand=operand+c;
             }
             return operand;

}

long double function_to_double(std::string f, std::vector<Variable> parameters)
{
    std::string expression, c;
    int i=0, limit=f.length();
    while(i<limit){//name
        if(is_operator_or_parenthesis(f[i]))
        {
            expression=expression+f[i];
            i++;
            continue;
        }
        read_operand(f, c, i);
        bool encontrou=false;

        for(Variable variable:parameters){
            if(variable.name==c)
            {
                std::stringstream ss;
                std::string aux;
                ss<<std::setprecision(400)<<variable.value;
                ss>>std::setprecision(400)>>aux;
                expression=expression+aux;
                encontrou=true;
                break;
            }
        }
        if(encontrou==false && c=="e")
        {
            expression=expression+euler;
        }
        else if(encontrou==false && c=="p")
        {
            expression=expression+pi;
        }
        else if(encontrou==false)
        {
            expression=expression+c;
        }

    }
return expression_evaluator(expression);
}


long double partial(std::string f, std::string variable, std::vector<Variable> parameters)//name
{
    const static long double h=0.00012207031;

    int i=0;
    for(; i<parameters.size(); ++i){
        if(parameters[i].name==variable)
        {
            break;
        }

    }
    long double auxiliar1, auxiliar2;
    auxiliar1=function_to_double(f, parameters);
    parameters[i].value+=h;
    auxiliar2=function_to_double(f, parameters);
    return (auxiliar2-auxiliar1)/h;

}

void read_parameters(std::fstream& file, std::vector<Variable>&parameters){

std::string name;//name
long double value,sigma;
while(file>>name>>value>>sigma)
{
    parameters.push_back({name, value, sigma});
}
return;
}

long double uncertainty(std::string f, std::vector<Variable> parameters){

long double var=0;
for(Variable parameter : parameters){

    var+=parameter.sigma*parameter.sigma*partial(f, parameter.name, parameters)*partial(f, parameter.name, parameters);
}
return std::sqrt(var);
}

void read_variable(std::fstream& file, std::vector<Variable>& parameters){


std::string name;//name
long double n=0;
long double equip_sigma, measure, sum=0, sum_squares=0, mean, sigma, standard_deviation=0, standard_deviation_mean=0;
file>>name>>equip_sigma;

while(file>>measure)
{
    if(measure==0)
    {
        break;
    }
    sum+=measure;
    sum_squares+=measure*measure;
    n+=1;
}
if(n==1)
{
    mean=sum;
    parameters.push_back({name, mean, equip_sigma});
    return;
}

mean=sum/n;
standard_deviation=std::sqrt((sum_squares-mean*mean*n)/(n-1));
standard_deviation_mean=std::sqrt((sum_squares-mean*mean*n)/(n*(n-1)));//already takes sigma of the mean
sigma=std::sqrt(standard_deviation_mean*standard_deviation_mean+equip_sigma*equip_sigma);//accounts for equip uncertainty
parameters.push_back({name, mean, sigma, standard_deviation, standard_deviation_mean});
return;
}

void read_txt_table(std::fstream& file, std::vector<Variable>& parameters){

while(file.peek()!='*'){
    read_variable(file, parameters);

}
return;
}

void read_function(std::fstream& file, std::vector<std::string>& functions){

char star;//we will discart the * character from the 0* sequence
std::string f;
file>>star;
while(file>>f){

functions.push_back(f);

}
return;
}

void click_anything(){

char c;
std::cout<<"\n\nclick any key to leave to leave\n";

std::cin>>std::noskipws>>c;
return;
}

void receptionist(std::string& filename, std::string& command){

std::cout<<"enter file name with table:\n\n";
std::cin>>filename;
if(std::getline(std::cin, command))
{

}
return;
}

std::string get_edited_table_name(std::string filename){

std::string new_name;
int len=filename.length();
for(int i=0; i<len; ++i)
{
    if(filename[i]=='.' &&filename[i+1]=='t' &&filename[i+2]=='x' &&filename[i+3]=='t')
    {
        break;
    }
    new_name=new_name+filename[i];
}
return new_name+"_final.txt";
}

void write(std::fstream& output, std::string command, std::vector<Variable> parameters, std::vector<std::string> functions){

output<<std::setprecision(40);

if(command.find("linear_fit") != std::string::npos)
{
    long double angular_value, angular_sigma, linear_value, linear_sigma, error_squared;
    angular_value=angular_right(parameters, parameters.size()/2);
    angular_sigma=angular_sigma_right(parameters, parameters.size()/2);
    linear_value=linear_right(parameters, parameters.size()/2);
    linear_sigma=linear_sigma_right(parameters, parameters.size()/2);
    error_squared=error(parameters);

    parameters.push_back({"angular", angular_value, angular_sigma});
    parameters.push_back({"linear", linear_value, linear_sigma});
    parameters.push_back({"error", error_squared});
}

if(command.find("quadratic_fit") != std::string::npos)
{
    long double a2, a2_sigma, a1, a1_sigma, a0, a0_sigma, quadratic_error_squared, n=parameters.size()/2;
    a2=a_2(parameters, n);
    a2_sigma=a_2_sigma(parameters, n);
    a1=a_1(parameters, n);
    a1_sigma=a_1_sigma(parameters, n);
    a0=a_0(parameters, n);
    a0_sigma=a_0_sigma(parameters, n);
    quadratic_error_squared=error_parabola(parameters);

    parameters.push_back({"a_2", a2, a2_sigma});
    parameters.push_back({"a_1", a1, a1_sigma});
    parameters.push_back({"a_0", a0, a0_sigma});
    parameters.push_back({"error", quadratic_error_squared});
}

for(Variable variavel:parameters){
    if(variavel.name=="angular" || variavel.name=="a_2")
    {
        output<<'\n';
    }

    output<<variavel.name<<' '<<variavel.value<<" ± "<<variavel.sigma<<'\n';

}

for(std::string f:functions){

output<<'\n'<<f<<"\n\n";
output<<function_to_double(f, parameters)<<" ± "<<uncertainty(f, parameters)<<'\n';

}

return;

}

long double angular_right(std::vector<Variable> parameters, long double n){
    long double sum_xy=0, sum_x=0, sum_y=0, sum_x2=0;

    for(int i=0; i<n; ++i){
        sum_xy+=parameters[i].value*parameters[i+n].value;
        sum_x+=parameters[i].value;
        sum_y+=parameters[i+n].value;
        sum_x2+=parameters[i].value*parameters[i].value;
    }

    long double angular=(n*sum_xy-sum_x*sum_y)/(n*sum_x2-sum_x*sum_x);
    return angular;
}

long double linear_right(std::vector<Variable> parameters, long double n){
    long double sum_xy=0, sum_x=0, sum_y=0, sum_x2=0;

    for(int i=0; i<n; ++i){
        sum_xy+=parameters[i].value*parameters[i+n].value;
        sum_x+=parameters[i].value;
        sum_y+=parameters[i+n].value;
        sum_x2+=parameters[i].value*parameters[i].value;
    }

    long double linear=(sum_x2*sum_y-sum_x*sum_xy)/(n*sum_x2-sum_x*sum_x);
    return linear;
}

long double angular_sigma_right(std::vector<Variable> parameters, long double n){
    const static long double h=0.00012207031;
    long double var=0,f1, f2;
    for(int i=0; i<2*n; ++i){
        f1=angular_right(parameters, n);
        parameters[i].value+=h;
        f2=angular_right(parameters, n);
        parameters[i].value-=h;
        var+=parameters[i].sigma*((f2-f1)/h)*parameters[i].sigma*((f2-f1)/h);
    }
    return std::sqrt(var);
}

long double linear_sigma_right(std::vector<Variable> parameters, long double n){
    const static long double h=0.00012207031;
    long double var=0,f1, f2;
    for(int i=0; i<2*n; ++i){
        f1=linear_right(parameters, n);
        parameters[i].value+=h;
        f2=linear_right(parameters, n);
        parameters[i].value-=h;
        var+=parameters[i].sigma*((f2-f1)/h)*parameters[i].sigma*((f2-f1)/h);
    }
    return std::sqrt(var);
}

long double error(std::vector<Variable> parameters){
    const int n=parameters.size()/2;
    long double X=0, a=angular_right(parameters, n), b=linear_right(parameters, n);
    for(int i=0; i<n; ++i){
        X+=(parameters[i+n].value-(a*parameters[i].value+b))*(parameters[i+n].value-(a*parameters[i].value+b))/*/(parameters[i].sigma*parameters[i].sigma+parameters[i+n].sigma*parameters[i+n].sigma)*/;
    }
    return X;
}

long double a_2(std::vector<Variable> parameters, long double n){
    long double s0=0, s1=0, s2=0, s3=0, s4=0, yx0=0, yx1=0, yx2=0;
    for(int i=0; i<n; ++i){
        s0+=1;
        s1+=parameters[i].value;
        s2+=std::pow(parameters[i].value, 2);
        s3+=std::pow(parameters[i].value, 3);
        s4+=std::pow(parameters[i].value, 4);
        yx0+=parameters[i+n].value;
        yx1+=parameters[n+i].value*parameters[i].value;
        yx2+=parameters[n+i].value*std::pow(parameters[i].value, 2);
    }
    return (yx0*s1*s1*s3-yx0*s1*s2*s2+yx1*s1*s1*s2-yx1*s0*s1*s3+yx2*s0*s1*s2-yx2*s1*s1*s1)/(s0*s1*s2*s4-s0*s1*s3*s3+s1*s1*s2*s3-s1*s1*s1*s4+s1*s1*s2*s3-s1*s2*s2*s2);
}

long double a_1(std::vector<Variable> parameters, long double n){
    long double s0=0, s1=0, s2=0, s3=0, s4=0, yx0=0, yx1=0, yx2=0;
    for(int i=0; i<n; ++i){
        s0+=1;
        s1+=parameters[i].value;
        s2+=std::pow(parameters[i].value, 2);
        s3+=std::pow(parameters[i].value, 3);
        s4+=std::pow(parameters[i].value, 4);
        yx0+=parameters[i+n].value;
        yx1+=parameters[n+i].value*parameters[i].value;
        yx2+=parameters[n+i].value*std::pow(parameters[i].value, 2);
    }
    return (yx2*s1-yx1*s2-a_2(parameters, n)*(s1*s4-s2*s3))/(s1*s3-s2*s2);
}
long double a_0(std::vector<Variable> parameters, long double n){
    long double s0=0, s1=0, s2=0, s3=0, s4=0, yx0=0, yx1=0, yx2=0;
    for(int i=0; i<n; ++i){
        s0+=1;
        s1+=parameters[i].value;
        s2+=std::pow(parameters[i].value, 2);
        s3+=std::pow(parameters[i].value, 3);
        s4+=std::pow(parameters[i].value, 4);
        yx0+=parameters[i+n].value;
        yx1+=parameters[n+i].value*parameters[i].value;
        yx2+=parameters[n+i].value*std::pow(parameters[i].value, 2);
    }
    return (yx0-s1*a_1(parameters, n)-s2*a_2(parameters, n))/(s0);
}

long double a_2_sigma(std::vector<Variable> parameters, long double n){
    const static long double h=0.00012207031;
    long double var=0, f2, f1;
    for(int i=0; i<2*n; i++){
        parameters[i].value+=h;
        f2=a_2(parameters, n);
        parameters[i].value-=h;
        f1=a_2(parameters, n);
        var+=((f2-f1)/h)*((f2-f1)/h)*parameters[i].sigma*parameters[i].sigma;
    }
    return std::sqrt(var);
}

long double a_1_sigma(std::vector<Variable> parameters, long double n){
    const static long double h=0.00012207031;
    long double var=0, f2, f1;
    for(int i=0; i<2*n; i++){
        parameters[i].value+=h;
        f2=a_1(parameters, n);
        parameters[i].value-=h;
        f1=a_1(parameters, n);
        var+=((f2-f1)/h)*((f2-f1)/h)*parameters[i].sigma*parameters[i].sigma;
    }
    return std::sqrt(var);
}

long double a_0_sigma(std::vector<Variable> parameters, long double n){
    const static long double h=0.00012207031;
    long double var=0, f2, f1;
    for(int i=0; i<2*n; i++){
        parameters[i].value+=h;
        f2=a_0(parameters, n);
        parameters[i].value-=h;
        f1=a_0(parameters, n);
        var+=((f2-f1)/h)*((f2-f1)/h)*parameters[i].sigma*parameters[i].sigma;
    }
    return std::sqrt(var);
}

long double error_parabola(std::vector<Variable> parameters){
    const int n=parameters.size()/2;
    long double X=0, a2=a_2(parameters, n), a1=a_1(parameters, n), a0=a_0(parameters, n);
    for(int i=0; i<n; ++i){
        X+=(parameters[i+n].value-(a2*parameters[i].value*parameters[i].value+a1*parameters[i].value+a0))*(parameters[i+n].value-(a2*parameters[i].value*parameters[i].value+a1*parameters[i].value+a0))/*/(parameters[i].sigma*parameters[i].sigma+parameters[i+n].sigma*parameters[i+n].sigma)*/;
    }
    return X;
}

bool is_operator_or_parenthesis(char c){

if(c=='+')
{
    return true;
}

if(c=='-')
{
    return true;
}

if(c=='*')
{
    return true;
}

if(c=='/')
{
    return true;
}

if(c=='^')
{
    return true;
}

if(c=='S')
{
    return true;
}

if(c=='C')
{
    return true;
}

if(c=='T')
{
    return true;
}

if(c=='L')
{
    return true;
}

if(c=='(')
{
    return true;
}

if(c==')')
{
    return true;
}

return false;

}

void read_operand(std::string& function, std::string& c, int& i){
c="";
int limit=function.length();
while(i<limit && !is_operator_or_parenthesis(function[i])){
c=c+function[i];
i++;
}
return;
}
