 #include "graph.h"
string reverse_path(string path)
{
	string reverse_path("");
	char buffer[BUFFER_LENTH];
	char dst[MAX_K][BUFFER_LENTH];
	strncpy(buffer, path.c_str(), path.length()+1);
	//buffer[path.length + 1] = '\0';
	int cnt = split(dst, buffer, " ");
	for (int i = cnt-1; i >= 0; i--)
	{
		reverse_path += (dst[i]);
		reverse_path += " ";
	}
	return reverse_path;
}


bool circle_detection_string(string path)// true means there is a circle
{
	map<string, bool> test;
	char buffer[BUFFER_LENTH];
	char dst[MAX_K][BUFFER_LENTH];
	strncpy(buffer, path.c_str(), path.length() + 1);
	//buffer[path.length + 1] = '\0';
	int cnt = split(dst, buffer, " ");
	for (int i = cnt - 1; i >= 0; i--)
	{
		if (test.find(dst[i]) == test.end())
		{
			test.insert(pair<string, bool>(dst[i], true));
		}
		else
		{
			return true;
		}
	}
	return false;
}

bool circle_detection(vector<NODE_TYPE> path)// true means there is a circle
{
	map<NODE_TYPE, bool> test;
	for (auto iter = path.begin(); iter != path.end(); iter++)
	{
		if (test.find(*iter) == test.end())
		{
			test.insert(pair<NODE_TYPE, bool>( *iter,true));
		}
		else
		{
			return true;
		}
	}
	return false;
}
int split(char sdt[][BUFFER_LENTH], char* str, const char* spl)
{
	int n = 0;
	char * result = NULL;
	result = strtok(str, spl);
	while (result != NULL)
	{
		strcpy(sdt[n++], result);
		result = strtok(NULL, spl);
	}
	return n;
}

time_t StringToDatetime(string str)
{
	char *cha = (char*)str.data();             // 将string转换成char*。
	tm tm_;                                    // 定义tm结构体。
	int year, month, day, hour, minute, second;// 定义时间的各个int临时变量。
	sscanf(cha, "%d-%d-%d %d:%d:%d", &year, &month, &day, &hour, &minute, &second);// 将string存储的日期时间，转换为int临时变量。
	tm_.tm_year = year - 1900;                 // 年，由于tm结构体存储的是从1900年开始的时间，所以tm_year为int临时变量减去1900。
	tm_.tm_mon = month - 1;                    // 月，由于tm结构体的月份存储范围为0-11，所以tm_mon为int临时变量减去1。
	tm_.tm_mday = day;                         // 日。
	tm_.tm_hour = hour;                        // 时。
	tm_.tm_min = minute;                       // 分。
	tm_.tm_sec = second;                       // 秒。
	tm_.tm_isdst = 0;                          // 非夏令时。
	time_t t_ = mktime(&tm_);                  // 将tm结构体转换成time_t格式。
	return t_;                                 // 返回值。 
}

string DatetimeToString(time_t time)
{
	tm *tm_ = localtime(&time);                // 将time_t格式转换为tm结构体
	int year, month, day, hour, minute, second;// 定义时间的各个int临时变量。
	year = tm_->tm_year + 1900;                // 临时变量，年，由于tm结构体存储的是从1900年开始的时间，所以临时变量int为tm_year加上1900。
	month = tm_->tm_mon + 1;                   // 临时变量，月，由于tm结构体的月份存储范围为0-11，所以临时变量int为tm_mon加上1。
	day = tm_->tm_mday;                        // 临时变量，日。
	hour = tm_->tm_hour;                       // 临时变量，时。
	minute = tm_->tm_min;                      // 临时变量，分。
	second = tm_->tm_sec;                      // 临时变量，秒。
	char yearStr[5], monthStr[3], dayStr[3], hourStr[3], minuteStr[3], secondStr[3];// 定义时间的各个char*变量。
	sprintf(yearStr, "%d", year);              // 年。
	sprintf(monthStr, "%d", month);            // 月。
	sprintf(dayStr, "%d", day);                // 日。
	sprintf(hourStr, "%d", hour);              // 时。
	sprintf(minuteStr, "%d", minute);          // 分。
	if (minuteStr[1] == '\0')                  // 如果分为一位，如5，则需要转换字符串为两位，如05。
	{
		minuteStr[2] = '\0';
		minuteStr[1] = minuteStr[0];
		minuteStr[0] = '0';
	}
	sprintf(secondStr, "%d", second);          // 秒。
	if (secondStr[1] == '\0')                  // 如果秒为一位，如5，则需要转换字符串为两位，如05。
	{
		secondStr[2] = '\0';
		secondStr[1] = secondStr[0];
		secondStr[0] = '0';
	}
	char s[20];                                // 定义总日期时间char*变量。
	sprintf(s, "%s-%s-%s %s:%s:%s", yearStr, monthStr, dayStr, hourStr, minuteStr, secondStr);// 将年月日时分秒合并。
	string str(s);                             // 定义string变量，并将总日期时间char*变量作为构造函数的参数传入。
	return str;                                // 返回转换日期时间后的string变量。
}

void data_clean()
{
	int table_range[6] = {0,7,11,11,11,7};//table from 1 to 5 and range from 0 to table range[i]
	int table = 1;
	NODE_TYPE x, y;
	string table_name;
	string out_table_name;
	map<string, int> reordered_edges;
	NODE_TYPE cur_id = 0;
	ifstream inEdges;
	ofstream outEdges;
	char buffer[BUFFER_LENTH];
	//string x,y;
	char dst[20][BUFFER_LENTH];
	for (table = 1; table <= 5; table++)
	{
		for (int cur_index = 0; cur_index <= table_range[table]; cur_index++)
		{

			table_name = "table"+ to_string(table) +"/table" + to_string(table)+ "_"+ to_string(cur_index) + ".txt";
			out_table_name = "table" + to_string(table) + "/table" + to_string(table) + "_" + to_string(cur_index) + "out.txt";
			// cout << table_name << " is our talbe name" << endl;
			//for each table, input all the key-value pairs. Store and output all ordered key-value pairs.
			inEdges.open(table_name);
			outEdges.open(out_table_name);
			if (!inEdges.is_open())
			{
				// cout << "Error opening file : "<< table_name << endl;
				continue;
			}
			while (!inEdges.eof())
			{
				inEdges.getline(buffer, BUFFER_LENTH);
				int cnt = split(dst, buffer, ",");
				if (reordered_edges.find(dst[1]) == reordered_edges.end()) {
					reordered_edges[dst[1]] = cur_id;
					x = cur_id;
					cur_id++;
				}
				else {
					x = reordered_edges[dst[1]];
				}
				if (reordered_edges.find(dst[3]) == reordered_edges.end()) {
					reordered_edges[dst[3]] = cur_id;
					y = cur_id;
					cur_id++;
				}
				else {
					y = reordered_edges[dst[3]];
				}
				outEdges << x << " " << y << endl;
			}
			inEdges.close();
			outEdges.close();
		}
	}
}

NODE_TYPE find_max_node()
{
	NODE_TYPE max_node = 0;
	int table_range[6] = { 0,7,11,11,11,7 };//table from 1 to 5 and range from 0 to table range[i]
	int table = 1;
	NODE_TYPE x, y;
	string table_name;
	ifstream inEdges;
	char buffer[BUFFER_LENTH];
	//string x,y;
	for (table = 1; table <= 5; table++)
	{
		for (int cur_index = 0; cur_index <= table_range[table]; cur_index++)
		{

			table_name = "table" + to_string(table) + "/table" + to_string(table) + "_" + to_string(cur_index) + "out.txt";
			// cout << table_name << " is our talbe name" << endl;
			//for each table, input all the key-value pairs. Store and output all ordered key-value pairs.
			inEdges.open(table_name);
			if (!inEdges.is_open())
			{
				// cout << "Error opening file : " << table_name << endl;
				continue;
			}
			while (!inEdges.eof())
			{
				inEdges.getline(buffer, BUFFER_LENTH);
				extractEdges(buffer, x, y);
				max_node = x > max_node ? x : max_node;
				max_node = y > max_node ? y : max_node;
			}
			inEdges.close();
		}
	}
	// cout << max_node << " is our max_node " << endl;
	return max_node;
}