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
	char *cha = (char*)str.data();             // ��stringת����char*��
	tm tm_;                                    // ����tm�ṹ�塣
	int year, month, day, hour, minute, second;// ����ʱ��ĸ���int��ʱ������
	sscanf(cha, "%d-%d-%d %d:%d:%d", &year, &month, &day, &hour, &minute, &second);// ��string�洢������ʱ�䣬ת��Ϊint��ʱ������
	tm_.tm_year = year - 1900;                 // �꣬����tm�ṹ��洢���Ǵ�1900�꿪ʼ��ʱ�䣬����tm_yearΪint��ʱ������ȥ1900��
	tm_.tm_mon = month - 1;                    // �£�����tm�ṹ����·ݴ洢��ΧΪ0-11������tm_monΪint��ʱ������ȥ1��
	tm_.tm_mday = day;                         // �ա�
	tm_.tm_hour = hour;                        // ʱ��
	tm_.tm_min = minute;                       // �֡�
	tm_.tm_sec = second;                       // �롣
	tm_.tm_isdst = 0;                          // ������ʱ��
	time_t t_ = mktime(&tm_);                  // ��tm�ṹ��ת����time_t��ʽ��
	return t_;                                 // ����ֵ�� 
}

string DatetimeToString(time_t time)
{
	tm *tm_ = localtime(&time);                // ��time_t��ʽת��Ϊtm�ṹ��
	int year, month, day, hour, minute, second;// ����ʱ��ĸ���int��ʱ������
	year = tm_->tm_year + 1900;                // ��ʱ�������꣬����tm�ṹ��洢���Ǵ�1900�꿪ʼ��ʱ�䣬������ʱ����intΪtm_year����1900��
	month = tm_->tm_mon + 1;                   // ��ʱ�������£�����tm�ṹ����·ݴ洢��ΧΪ0-11��������ʱ����intΪtm_mon����1��
	day = tm_->tm_mday;                        // ��ʱ�������ա�
	hour = tm_->tm_hour;                       // ��ʱ������ʱ��
	minute = tm_->tm_min;                      // ��ʱ�������֡�
	second = tm_->tm_sec;                      // ��ʱ�������롣
	char yearStr[5], monthStr[3], dayStr[3], hourStr[3], minuteStr[3], secondStr[3];// ����ʱ��ĸ���char*������
	sprintf(yearStr, "%d", year);              // �ꡣ
	sprintf(monthStr, "%d", month);            // �¡�
	sprintf(dayStr, "%d", day);                // �ա�
	sprintf(hourStr, "%d", hour);              // ʱ��
	sprintf(minuteStr, "%d", minute);          // �֡�
	if (minuteStr[1] == '\0')                  // �����Ϊһλ����5������Ҫת���ַ���Ϊ��λ����05��
	{
		minuteStr[2] = '\0';
		minuteStr[1] = minuteStr[0];
		minuteStr[0] = '0';
	}
	sprintf(secondStr, "%d", second);          // �롣
	if (secondStr[1] == '\0')                  // �����Ϊһλ����5������Ҫת���ַ���Ϊ��λ����05��
	{
		secondStr[2] = '\0';
		secondStr[1] = secondStr[0];
		secondStr[0] = '0';
	}
	char s[20];                                // ����������ʱ��char*������
	sprintf(s, "%s-%s-%s %s:%s:%s", yearStr, monthStr, dayStr, hourStr, minuteStr, secondStr);// ��������ʱ����ϲ���
	string str(s);                             // ����string����������������ʱ��char*������Ϊ���캯���Ĳ������롣
	return str;                                // ����ת������ʱ����string������
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