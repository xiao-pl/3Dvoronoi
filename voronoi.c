#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>


// 存储点
typedef struct GroupPoint {
    uint32_t index;
    int point[3];             // 存储整型坐标点
    uint32_t neighbor;      // 邻居信息
    uint16_t region;        // 所属区域编号
    struct GroupPoint* next;
} GroupPoint;

// ========================
// queue
// ========================

// 结构体表示队列节点
typedef struct QueueNode {
    GroupPoint* data;
    struct QueueNode* next;
} QueueNode;

// 结构体表示队列
typedef struct {
    QueueNode* front;
    QueueNode* rear;
} Queue;

// 初始化队列
Queue* initQueue()
{
    Queue* queue = (Queue*)malloc(sizeof(Queue));
    if(queue != NULL) {
        queue->front = NULL;
        queue->rear = NULL;
    }
    return queue;
}

// 入队操作
void enqueue(Queue* queue, GroupPoint* point)
{
    QueueNode* newNode = (QueueNode*)malloc(sizeof(QueueNode));
    if(newNode != NULL) {
        newNode->data = point;
        newNode->next = NULL;

        if(queue->rear == NULL) {
            queue->front = newNode;
            queue->rear = newNode;
        }
        else {
            queue->rear->next = newNode;
            queue->rear = newNode;
        }
    }
}

// 出队操作
GroupPoint* dequeue(Queue* queue)
{
    if(queue->front == NULL) {
        return NULL; // 队列为空
    }

    QueueNode* temp = queue->front;
    GroupPoint* data = temp->data;

    queue->front = queue->front->next;
    if(queue->front == NULL) {
        queue->rear = NULL;
    }

    free(temp);
    return data;
}

// ========================
// dict
// ========================

// 结构体表示哈希表中的节点
typedef struct HashNode {
    GroupPoint* key;
    GroupPoint* value[];
} HashNode;

// 向哈希表中插入元素
void insert_value_to_hash_table(HashNode* hash_table[], GroupPoint* const key, GroupPoint* const value, uint32_t arr_point_len)
{
    unsigned int index = key->index;

    // 插入到哈希表中
    if(hash_table[index] == NULL) {
        // 创建新节点
        HashNode* new_node = (HashNode*)malloc(sizeof(HashNode) + arr_point_len * sizeof(GroupPoint*));
        for(uint32_t i = 0; i < arr_point_len; ++i) {
            new_node->value[i] = NULL;
        }
        new_node->key = key;
        new_node->value[value->index] = value;
        hash_table[index] = new_node;
    }
    else {
        HashNode* new_node = hash_table[index];
        new_node->value[value->index] = value;
    }
}

// 向哈希表中插入元素
void insert_key_to_hash_table(HashNode* hash_table[], GroupPoint* const key, uint32_t arr_point_len)
{
    unsigned int index = key->index;

    // 插入到哈希表中
    if(hash_table[index] == NULL) {
        // 创建新节点
        HashNode* new_node = (HashNode*)malloc(sizeof(HashNode) + arr_point_len * sizeof(GroupPoint*));
        for(uint32_t i = 0; i < arr_point_len; ++i) {
            new_node->value[i] = NULL;
        }
        new_node->key = key;
        hash_table[index] = new_node;
    }
}



// ========================
// hash
// ========================


#define HASH_TABLE_SIZE 10000 // 哈希表大小

// 3D点的结构体
typedef struct {
    float x;
    float y;
    float z;
} Point3D;

// 哈希节点的结构体
typedef struct HashNode_point {
    Point3D point;
    struct HashNode_point* next;
} HashNode_point;

// 哈希表的结构体
typedef struct {
    HashNode_point* buckets[HASH_TABLE_SIZE];
} HashTable;

// 计算点的哈希值
uint32_t hash_function(float x, float y, float z)
{
    // 这里可以根据具体情况设计哈希函数
    return (uint32_t)(x * 1000 + y * 100 + z);
}

// 初始化哈希表
void init_hash_table(HashTable* table)
{
    for(int i = 0; i < HASH_TABLE_SIZE; ++i) {
        table->buckets[i] = NULL;
    }
}

// 向哈希表中插入点云数据
void insert_point(HashTable* table, Point3D point)
{
    uint32_t index = hash_function(point.x, point.y, point.z) % HASH_TABLE_SIZE;

    HashNode_point* newNode = (HashNode_point*)malloc(sizeof(HashNode_point));
    newNode->point = point;
    newNode->next = table->buckets[index];
    table->buckets[index] = newNode;
}

// 判断点是否在点云中
int is_point_in_cloud(HashTable* table, float x, float y, float z)
{
    uint32_t index = hash_function(x, y, z) % HASH_TABLE_SIZE;

    HashNode_point* current = table->buckets[index];
    while(current != NULL) {
        if(current->point.x == x && current->point.y == y && current->point.z == z) {
            return 1; // 找到了点
        }
        current = current->next;
    }
    return 0; // 没有找到点
}

// =================================================================================

GroupPoint m_groupPoint;    // 所有点
uint32_t number_of_voxels = 0;

// =================================================================================

/**
 * @brief 拷贝点
 *
 * @param src_point 源点
 * @param dst_point 目标点
 */
void copy_point(GroupPoint* src_point, GroupPoint* dst_point)
{
    dst_point->index = src_point->index;
    for(int i = 0;i < 3;++i) {
        dst_point->point[i] = src_point->point[i];
    }
    dst_point->neighbor = src_point->neighbor;
    dst_point->next = src_point->next;
    dst_point->region = src_point->region;
}

// 打印进度条的函数
void print_progress(float progress)
{
    int bar_length = 50; // 进度条长度
    int current_length = progress * bar_length; // 根据进度计算当前长度

    printf("\rProgress: [");
    for(int i = 0; i < bar_length; ++i) {
        if(i < current_length) {
            printf("#");
        }
        else {
            printf(" ");
        }
    }
    printf("] %.2f%%", progress * 100); // 打印百分比

    fflush(stdout); // 立即刷新输出缓冲区
}

/**
 * @brief 获取GroupPoint的长度
 *
 * @param points
 * @return uint32_t
 */
uint32_t get_GroupPoint_len(GroupPoint* points)
{
    uint32_t len = 0;
    while(points) {
        ++len;
        points = points->next;
    }
    return len;
}

/**
 * @brief 读取点云数据
 *
 * @param filename *.txt点云文件路径
 * @param points 所有点
 */
int parse_input(const char* filename, GroupPoint* const points)
{
    FILE* file;
    char buffer[1000];

    GroupPoint* current_groupPoint = points;
    GroupPoint* last_groupPoint = NULL;

    file = fopen(filename, "r");
    if(file == NULL) {
        printf("Unable to open the file.\n");
        return 1;
    }

    number_of_voxels = 0;

    // 读取文件内容
    while(fgets(buffer, sizeof(buffer), file) != NULL) {
        // 每行三个空格隔开的整数值(x y z)
        int x, y, z;
        sscanf_s(buffer, "%d %d %d", &x, &y, &z);

        current_groupPoint->index = number_of_voxels++;
        current_groupPoint->point[0] = x;
        current_groupPoint->point[1] = y;
        current_groupPoint->point[2] = z;
        current_groupPoint->neighbor = 0;
        current_groupPoint->region = 0;

        GroupPoint* newNode = (GroupPoint*)malloc(sizeof(GroupPoint));
        current_groupPoint->next = newNode;
        last_groupPoint = current_groupPoint;
        current_groupPoint = newNode;
    }

    if(last_groupPoint) {
        free(current_groupPoint);
        last_groupPoint->next = NULL;
    }

    fclose(file);

    return 0;
}

/**
 * @brief 判断是否在点组中
 *
 * @param points
 * @param x
 * @param y
 * @param z
 * @return int
 */
int is_in_GroupPoint(GroupPoint* points, int x, int y, int z)
{
    while(points) {
        if(points->point[0] == x && points->point[1] == y && points->point[2] == z) {
            return 1;
        }

        points = points->next;
    }

    return 0;
}

/**
 * @brief 构造所有点的的邻居
 *
 * @param m_points 所有点
 * @param hash_table 检索表
 */
void mark_neighbors(GroupPoint* const m_points, HashTable* hash_table)
{
    GroupPoint* points = m_points;
    int process = 0;

    while(points) {
        int starting_point_x = points->point[0] - 1;
        int starting_point_y = points->point[1] - 1;
        int starting_point_z = points->point[2] - 1;

        int neighbour_no = 0;
        for(int i = 0; i < 3; ++i) {
            for(int j = 0; j < 3; ++j) {
                for(int k = 0; k < 3; ++k) {
                    int point_to_check_x = starting_point_x + i;
                    int point_to_check_y = starting_point_y + j;
                    int point_to_check_z = starting_point_z + k;

                    // if(is_in_GroupPoint(m_points, point_to_check_x, point_to_check_y, point_to_check_z) &&
                    if(is_point_in_cloud(hash_table, point_to_check_x, point_to_check_y, point_to_check_z) &&
                        (point_to_check_x != points->point[0] || point_to_check_y != points->point[1] || point_to_check_z != points->point[2])) {          // 寻找九宫格中的邻居点
                        points->neighbor = points->neighbor | (1 << neighbour_no);
                    }

                    ++neighbour_no;
                }
            }
        }

        points = points->next;
        ++process;
        if(process % 20 == 0 || process == number_of_voxels) {
            print_progress(process * 1.0 / number_of_voxels);
        }
    }
    printf("\n");

}

/**
 * @brief 随机挑选种子点
 *
 * @param points 所有点
 * @param percent 挑选概率(0 ~ 100)
 * @return GroupPoint*
 */
GroupPoint* random_seeds(GroupPoint* points, double percent)
{
    GroupPoint* seed_points = NULL;

    while(points) {
        double random_number = (double)rand() / RAND_MAX * 100.0; // 生成 0 到 100 之间的随机浮点数

        // 根据随机数比例选择点
        if(random_number < percent) {
            // 添加符合条件的点到种子点集合中
            GroupPoint* newSeedPoint = (GroupPoint*)malloc(sizeof(GroupPoint));

            copy_point(points, newSeedPoint);

            newSeedPoint->next = seed_points;
            seed_points = newSeedPoint;
        }

        points = points->next;
    }

    return seed_points; // 返回种子点集合
}

/**
 * @brief 根据区域编号获取种子点
 *
 * @param region_number 区域编号
 * @param seed_points 种子点组
 * @return GroupPoint*
 */
GroupPoint* get_seed_point(int region_number, GroupPoint* seed_points)
{
    while(seed_points) {
        if(seed_points->region == region_number) {
            return seed_points;
        }

        seed_points = seed_points->next;
    }

    return NULL;
}

/**
 * @brief 根据索引获取点
 *
 * @param index 索引
 * @param points 所有点
 * @return GroupPoint*
 */
GroupPoint* get_point_by_index(GroupPoint* points, int index)
{
    while(points) {
        if(points->index == index) {
            return points;
        }
        points = points->next;
    }
    return NULL;
}

GroupPoint* get_point_by_xyz(GroupPoint* points, int x, int y, int z)
{
    while(points) {
        if(points->point[0] == x && points->point[1] == y && points->point[2] == z) {
            return points;
        }

        points = points->next;
    }

    return NULL;
}

/**
 * @brief 增量式BFS检索
 *
 * @param points 所有点
 * @param seed_points 种子点
 * @param adjacent_regions 种子点相邻的所有种子点, ！！！size应该为number_of_voxels！！！
 * @param boundary_points 种子点的所有边界点, 靠近分界线的点, ！！！size应该为number_of_voxels！！！
 * @return incremental_bfs_resulkt
 */
void incremental_bfs(GroupPoint* const points, GroupPoint* const seed_points, HashNode* adjacent_regions[], HashNode* boundary_points[])
{
    GroupPoint* t_points = points;
    while(t_points) {
        t_points->region = 0;       // 区域编号重置
        t_points = t_points->next;
    }

    int region_number = 1;          // 区域编号
    Queue* queue = initQueue();     // BFS查询队列

    GroupPoint* t_seed_points = seed_points;
    while(t_seed_points) {           // 每个种子点赋予唯一编号
        t_seed_points->region = region_number;
        t_points = get_point_by_index(points, t_seed_points->index);
        t_points->region = region_number;
        enqueue(queue, t_points);
        ++region_number;
        t_seed_points = t_seed_points->next;
    }

    GroupPoint* point;
    while((point = dequeue(queue)) != NULL) {
        region_number = point->region;

        int starting_point_x = point->point[0] - 1;
        int starting_point_y = point->point[1] - 1;
        int starting_point_z = point->point[2] - 1;

        int neighbour_no = -1;
        for(int i = 0; i < 3; ++i) {
            for(int j = 0; j < 3; ++j) {
                for(int k = 0; k < 3; ++k) {
                    ++neighbour_no;

                    if((1 << neighbour_no) & point->neighbor) {
                        int neighbour_x = starting_point_x + i;
                        int neighbour_y = starting_point_y + j;
                        int neighbour_z = starting_point_z + k;

                        GroupPoint* neighbour = get_point_by_xyz(points, neighbour_x, neighbour_y, neighbour_z);
                        if(neighbour->region != 0) {
                            if(neighbour->region != region_number) {     // 已被分配编号, 说明为边界点
                                GroupPoint* my_region = get_seed_point(region_number, seed_points);             // 获取当前区域的种子点
                                GroupPoint* adjacent_region = get_seed_point(neighbour->region, seed_points);   // 获取该边界点的种子点
                                insert_value_to_hash_table(boundary_points, my_region, point, number_of_voxels);
                                insert_value_to_hash_table(adjacent_regions, my_region, adjacent_region, number_of_voxels);
                            }
                            continue;
                        }
                        else {
                            neighbour->region = region_number;
                            enqueue(queue, neighbour);
                        }
                    }
                    else {
                        continue;
                    }
                }
            }
        }
    }

    free(queue);
}

/**
 * @brief 根据区域编号分类
 *
 * @param points 所有点
 * @param seed_points 种子点
 * @param region_classified_points 分类后的点组, ！！！size应该为number_of_voxels！！！
 */
void classify(GroupPoint* points, GroupPoint* const seed_points, HashNode* region_classified_points[])
{
    GroupPoint* t_seed_points = seed_points;
    while(t_seed_points) {
        insert_key_to_hash_table(region_classified_points, t_seed_points, number_of_voxels);
        t_seed_points = t_seed_points->next;
    }
    while(points) {
        GroupPoint* my_region = get_seed_point(points->region, seed_points);
        if(my_region != NULL) {
            insert_value_to_hash_table(region_classified_points, my_region, points, number_of_voxels);
        }
        points = points->next;
    }
}

/**
 * @brief 为种子点分配颜色
 *
 * @param seed_points 种子点
 * @param adjacent_regions 种子点相邻的所有种子点
 * @param region_colors 种子点分配的颜色值, ！！！size应该为number_of_voxels！！！
 */
void assign_colors(GroupPoint* seed_points, HashNode* adjacent_regions[], int region_colors[])
{
    for(int i = 0;i < number_of_voxels;++i) {
        region_colors[i] = 0;
    }

    while(seed_points) {
        HashNode* adjacent_region = adjacent_regions[seed_points->index];
        if(adjacent_region != NULL) {
            GroupPoint** neighbours = adjacent_region->value;
            for(int color = 1;color < 28; ++color) {
                int flag = 1;
                for(int i = 0; i < number_of_voxels; ++i) {
                    if(neighbours[i] != NULL) {
                        GroupPoint* neighbour = neighbours[i];
                        if(region_colors[neighbour->index] == color) {
                            flag = 0;
                            break;
                        }
                        else {
                            continue;
                        }
                    }
                }
                if(flag) {
                    region_colors[seed_points->index] = color;
                    break;
                }
            }
        }

        seed_points = seed_points->next;
    }
}

/**
 * @brief 向文件中写入obj数据
 *
 * @param color
 * @param point_x
 * @param point_y
 * @param point_z
 * @param count
 * @param f
 */
void print_cube(int color, float point_x, float point_y, float point_z, int count, FILE* f)
{
    float point1_x = point_x - 0.5, point1_y = point_y - 0.5, point1_z = point_z + 0.5;
    float point2_x = point_x - 0.5, point2_y = point_y - 0.5, point2_z = point_z - 0.5;
    float point3_x = point_x + 0.5, point3_y = point_y - 0.5, point3_z = point_z - 0.5;
    float point4_x = point_x + 0.5, point4_y = point_y - 0.5, point4_z = point_z + 0.5;
    float point5_x = point_x - 0.5, point5_y = point_y + 0.5, point5_z = point_z + 0.5;
    float point6_x = point_x + 0.5, point6_y = point_y + 0.5, point6_z = point_z + 0.5;
    float point7_x = point_x + 0.5, point7_y = point_y + 0.5, point7_z = point_z - 0.5;
    float point8_x = point_x - 0.5, point8_y = point_y + 0.5, point8_z = point_z - 0.5;

    fprintf(f, "v %.1f %.1f %.1f\n", point1_x, point1_y, point1_z);
    fprintf(f, "v %.1f %.1f %.1f\n", point2_x, point2_y, point2_z);
    fprintf(f, "v %.1f %.1f %.1f\n", point3_x, point3_y, point3_z);
    fprintf(f, "v %.1f %.1f %.1f\n", point4_x, point4_y, point4_z);
    fprintf(f, "v %.1f %.1f %.1f\n", point5_x, point5_y, point5_z);
    fprintf(f, "v %.1f %.1f %.1f\n", point6_x, point6_y, point6_z);
    fprintf(f, "v %.1f %.1f %.1f\n", point7_x, point7_y, point7_z);
    fprintf(f, "v %.1f %.1f %.1f\n", point8_x, point8_y, point8_z);
    fprintf(f, "usemtl color%d\n", color);
    fprintf(f, "f %d %d %d %d\n", count + 0, count + 1, count + 2, count + 3);
    fprintf(f, "f %d %d %d %d\n", count + 0, count + 3, count + 5, count + 4);
    fprintf(f, "f %d %d %d %d\n", count + 2, count + 6, count + 5, count + 3);
    fprintf(f, "f %d %d %d %d\n", count + 4, count + 5, count + 6, count + 7);
    fprintf(f, "f %d %d %d %d\n", count + 0, count + 4, count + 7, count + 1);
    fprintf(f, "f %d %d %d %d\n", count + 7, count + 6, count + 2, count + 1);
}


/**
 * @brief 生成可视化文件
 *
 * @param f 文件
 * @param seed_points 种子点
 * @param voronoi_cells 按区域编号分组后的点
 * @param region_colors 每个区域的颜色
 */
void generate_obj_file(FILE* f, GroupPoint* seed_points, HashNode* voronoi_cells[], int region_colors[])
{
    int count = 1;
    GroupPoint* t_seed_points = seed_points;
    while(t_seed_points) {
        if(region_colors[t_seed_points->index] != 0 && voronoi_cells[t_seed_points->index] != NULL) {
            GroupPoint** list_of_points = voronoi_cells[t_seed_points->index]->value;
            int region_color = region_colors[t_seed_points->index];
            for(int i = 0; i < number_of_voxels; ++i) {
                if(list_of_points[i] != NULL) {
                    GroupPoint* point = list_of_points[i];
                    if(is_in_GroupPoint(seed_points, point->point[0], point->point[1], point->point[2])) {
                        print_cube(0, point->point[0], point->point[1], point->point[2], count, f);
                    }
                    else {
                        print_cube(region_color, point->point[0], point->point[1], point->point[2], count, f);
                    }

                    count += 8;
                }
            }
        }

        t_seed_points = t_seed_points->next;
    }
}


float manhattan_distance(float p1_x, float p1_y, float p1_z, float p2_x, float p2_y, float p2_z)
{
    return fabs(p1_x - p2_x) + fabs(p1_y - p2_y) + fabs(p1_z - p2_z);
}

float get_max_distance_difference(int point_x, int point_y, int point_z, GroupPoint* boundary[])
{
    float min_distance = FLT_MAX;
    float max_distance = -FLT_MAX;

    for(int i = 0; i < number_of_voxels; ++i) {
        if(boundary[i] != NULL) {
            float distance = manhattan_distance(point_x, point_y, point_z, boundary[i]->point[0], boundary[i]->point[1], boundary[i]->point[2]);
            if(distance > max_distance) {
                max_distance = distance;
            }
            if(distance < min_distance) {
                min_distance = distance;
            }
        }
    }

    return max_distance - min_distance;
}

/**
 * @brief 生成新的种子点
 *
 * @param seed_points
 * @param voronoi_cells
 * @param boundary_points
 * @param new_seed_points 新的种子点, ！！！size为number_of_voxels, 且初始为NULL！！！
 */
void generate_new_seeds(GroupPoint* seed_points, HashNode* voronoi_cells[], HashNode* boundary_points[], GroupPoint* new_seed_points[])
{
    while(seed_points) {
        GroupPoint** region = voronoi_cells[seed_points->index]->value;
        GroupPoint** boundary = boundary_points[seed_points->index]->value;

        GroupPoint* new_seed = (GroupPoint*)malloc(sizeof(GroupPoint));
        copy_point(seed_points, new_seed);

        float min_difference = FLT_MAX;

        for(int i = 0; i < number_of_voxels; ++i) {
            if(region[i] != NULL) {
                GroupPoint* point = region[i];
                float difference = get_max_distance_difference(point->point[0], point->point[1], point->point[2], boundary);
                if(difference < min_difference) {           // 寻找在区域内更加居中的点作为新种子点
                    copy_point(point, new_seed);
                }
            }
        }

        GroupPoint* key_seed = (GroupPoint*)malloc(sizeof(GroupPoint));
        copy_point(seed_points, key_seed);
        new_seed_points[seed_points->index] = key_seed;
        new_seed->next = NULL;
        key_seed->next = new_seed;

        seed_points = seed_points->next;
    }
}

/**
 * @brief 改变种子点
 *
 * @param new_seed_points
 * @param region_colors
 * @param new_region_colors 新的区域色彩, ！！！size为number_of_voxels！！！
 * @return GroupPoint*
 */
GroupPoint* change_seeds(GroupPoint* new_seed_points[], int region_colors[], int new_region_colors[])
{
    GroupPoint* new_seeds = NULL;

    for(int i = 0; i < number_of_voxels; ++i) {
        if(new_seed_points[i] != NULL) {
            GroupPoint* seed = new_seed_points[i];
            GroupPoint* new_seed = seed->next;
            int color = region_colors[seed->index];
            new_region_colors[new_seed->index] = color;

            GroupPoint* newSeedPoint = (GroupPoint*)malloc(sizeof(GroupPoint));
            copy_point(new_seed, newSeedPoint);
            newSeedPoint->next = new_seeds;
            new_seeds = newSeedPoint;
        }
    }

    return new_seeds;
}

int find_color(int colors[], GroupPoint* neighbors[])
{
    int color = 1;

    for(int ii = 0; ii < 27; ++ii) {
        int flag = 0;
        for(int i = 0; i < number_of_voxels; ++i) {
            if(neighbors[i] != NULL) {
                GroupPoint* neighbor = neighbors[i];
                int neighbor_color = colors[neighbor->index];
                if(color == neighbor_color) {
                    flag = 1;
                    break;
                }
            }
        }

        if(flag)
            color += 1;
        else
            return color;
    }

    return color;
}

/**
 * @brief
 *
 * @param seed_points
 * @param region_colors
 * @param adjacent_regions
 * @param new_colors 新的区域色彩, ！！！size为number_of_voxels！！！
 */
void recolor(GroupPoint* seed_points, int region_colors[], HashNode* adjacent_regions[], int new_colors[])
{
    while(seed_points) {
        GroupPoint** neighbors = adjacent_regions[seed_points->index]->value;
        int my_color = region_colors[seed_points->index];
        int flag = 0;

        for(int i = 0;i < number_of_voxels;++i) {
            if(neighbors[i] != NULL) {
                GroupPoint* neighbor = neighbors[i];
                if(region_colors[neighbor->index] == my_color) {
                    flag = 1;
                    break;
                }
            }
        }
        if(flag) {
            my_color = find_color(region_colors, neighbors);
            region_colors[seed_points->index] = my_color;
        }
        new_colors[seed_points->index] = my_color;

        seed_points = seed_points->next;
    }
}

void free_GroupPoint(GroupPoint* f_point)
{
    while(f_point) {
        GroupPoint* t_point = f_point;
        f_point = f_point->next;
        free(t_point);
    }
}

void free_HashNode(HashNode* f_hashNode[])
{
    for(int i = 0; i < number_of_voxels; ++i) {
        if(f_hashNode[i] != NULL) {
            free(f_hashNode[i]);
        }
    }
    free(f_hashNode);
}

/**
 * @brief 最近点搜索函数
 * 
 * @param points 
 * @param numPoints 
 * @param queryPoint 
 * @return Point3D 
 */
Point3D findNearestPoint(GroupPoint* points, GroupPoint* seed_points, Point3D queryPoint)
{
    // 寻找给定点所属的区域
    Point3D ret = {0, 0, 0};

    int nearestRegionIndex = -1;
    float minDistance = INFINITY;

    while(points){
        float dist = manhattan_distance(points->point[0], points->point[1], points->point[2], queryPoint.x, queryPoint.y, queryPoint.z);
        if (dist < minDistance) {
            minDistance = dist;
            nearestRegionIndex = points->region;
        }
        points=points->next;
    }

    // 在找到所属区域后，寻找该区域内最近的点
    GroupPoint* my_region = get_seed_point(nearestRegionIndex, seed_points);
    ret.x = my_region->point[0];
    ret.y = my_region->point[1];
    ret.z = my_region->point[2];

    return ret;
}


// =================================================================================

/**
 * @brief 可视化测试
 *
 */
void test_0()
{
    char* input_filename = "D:/Desktop/voronoi/bunny.txt";
    int max_iterations = 2;

    printf("Reading input file\n");
    parse_input(input_filename, &m_groupPoint);
    printf("number_of_voxels: %d\n", number_of_voxels);
    if(number_of_voxels == 0) {
        return;
    }

    HashTable hash_table;
    init_hash_table(&hash_table);
    GroupPoint* t_point = &m_groupPoint;
    Point3D point_cloud;
    while(t_point) {
        point_cloud.x = t_point->point[0];
        point_cloud.y = t_point->point[1];
        point_cloud.z = t_point->point[2];
        insert_point(&hash_table, point_cloud);
        t_point = t_point->next;
    }


    printf("Marking neighbors\n");
    mark_neighbors(&m_groupPoint, &hash_table);

    printf("Initial random seeding\n");
    float seeds_num = 200;          // 设定种子点数量
    float seeds_percent = 100.0 * seeds_num / number_of_voxels;
    GroupPoint* seed_points = random_seeds(&m_groupPoint, seeds_percent);
    printf("seed points len: %d\n", get_GroupPoint_len(seed_points));

    printf("Starting iterations\n");
    int* region_colors = (int*)malloc(number_of_voxels * sizeof(int));
    memset(region_colors, 0, number_of_voxels * sizeof(int));
    for(int i = 0; i < max_iterations; ++i) {
        printf("Iteration number: %d\n", i + 1);

        HashNode** adjacent_regions = (HashNode**)malloc(number_of_voxels * sizeof(HashNode*));
        memset(adjacent_regions, 0, number_of_voxels * sizeof(HashNode*));
        HashNode** boundary_points = (HashNode**)malloc(number_of_voxels * sizeof(HashNode*));
        memset(boundary_points, 0, number_of_voxels * sizeof(HashNode*));
        printf("Incremental BFS\n");
        incremental_bfs(&m_groupPoint, seed_points, adjacent_regions, boundary_points);

        {
            int null_point = 0;
            GroupPoint* t_point = &m_groupPoint;
            while(t_point) {
                if(t_point->region == 0) {
                    ++null_point;
                }
                t_point = t_point->next;
            }
            printf("null point: %d\n", null_point);

        }

        HashNode** voronoi_cells = (HashNode**)malloc(number_of_voxels * sizeof(HashNode*));
        memset(voronoi_cells, 0, number_of_voxels * sizeof(HashNode*));
        printf("classify\n");
        classify(&m_groupPoint, seed_points, voronoi_cells);


        if(i == 0) {
            assign_colors(seed_points, adjacent_regions, region_colors);
        }
        else {
            int* new_colors = (int*)malloc(number_of_voxels * sizeof(int));
            memset(new_colors, 0, number_of_voxels * sizeof(int));
            recolor(seed_points, region_colors, adjacent_regions, new_colors);
            memcpy(region_colors, new_colors, number_of_voxels * sizeof(int));
            free(new_colors);
        }

        free_HashNode(adjacent_regions);

        printf("Create obj file\n");
        char out_file_name[100];
        strcpy(out_file_name, input_filename);
        out_file_name[strlen(input_filename) - 4] = '\0';
        sprintf(out_file_name, "%s%d.obj", out_file_name, i + 1);
        FILE* f = fopen(out_file_name, "w");
        fprintf(f, "mtllib ./colorfile.mtl\n");
        generate_obj_file(f, seed_points, voronoi_cells, region_colors);
        fclose(f);

        printf("Generating new seeds\n");
        GroupPoint** new_seed_points = (GroupPoint**)malloc(number_of_voxels * sizeof(GroupPoint*));
        memset(new_seed_points, 0, number_of_voxels * sizeof(GroupPoint*));
        generate_new_seeds(seed_points, voronoi_cells, boundary_points, new_seed_points);

        printf("Change seeds\n");
        int* new_colors = (int*)malloc(number_of_voxels * sizeof(int));
        memset(new_colors, 0, number_of_voxels * sizeof(int));
        GroupPoint* t_seed_points = change_seeds(new_seed_points, region_colors, new_colors);

        free_GroupPoint(seed_points);
        seed_points = t_seed_points;
        free(region_colors);
        region_colors = new_colors;

        free_HashNode(boundary_points);
        free_HashNode(voronoi_cells);
    }

    // *********** free *************
    free(region_colors);
    free_GroupPoint(m_groupPoint.next);
    free_GroupPoint(seed_points);
}



/**
 * @brief 最近邻搜索测试
 *
 */
void test_1()
{
    char* input_filename = "D:/Desktop/voronoi/bunny.txt";

    printf("Reading input file\n");
    parse_input(input_filename, &m_groupPoint);
    printf("number_of_voxels: %d\n", number_of_voxels);
    if(number_of_voxels == 0) {
        return;
    }

    HashTable hash_table;
    init_hash_table(&hash_table);
    GroupPoint* t_point = &m_groupPoint;
    Point3D point_cloud;
    while(t_point) {
        point_cloud.x = t_point->point[0];
        point_cloud.y = t_point->point[1];
        point_cloud.z = t_point->point[2];
        insert_point(&hash_table, point_cloud);
        t_point = t_point->next;
    }


    printf("Marking neighbors\n");
    mark_neighbors(&m_groupPoint, &hash_table);

    printf("Initial random seeding\n");
    float seeds_num = 200;          // 设定种子点数量
    float seeds_percent = 100.0 * seeds_num / number_of_voxels;
    GroupPoint* seed_points = random_seeds(&m_groupPoint, seeds_percent);
    printf("seed points len: %d\n", get_GroupPoint_len(seed_points));


    HashNode** adjacent_regions = (HashNode**)malloc(number_of_voxels * sizeof(HashNode*));
    memset(adjacent_regions, 0, number_of_voxels * sizeof(HashNode*));
    HashNode** boundary_points = (HashNode**)malloc(number_of_voxels * sizeof(HashNode*));
    memset(boundary_points, 0, number_of_voxels * sizeof(HashNode*));
    printf("Incremental BFS\n");
    incremental_bfs(&m_groupPoint, seed_points, adjacent_regions, boundary_points);


    HashNode** voronoi_cells = (HashNode**)malloc(number_of_voxels * sizeof(HashNode*));
    memset(voronoi_cells, 0, number_of_voxels * sizeof(HashNode*));
    printf("classify\n");
    classify(&m_groupPoint, seed_points, voronoi_cells);


    // **************** 最近邻搜索测试 ********************//
    Point3D queryPoint = {-56,   80  , 15};         // 测试点
    Point3D nearestPoint = findNearestPoint(&m_groupPoint, seed_points, queryPoint);

    // 最近点的坐标
    printf("Nearest point: (%f, %f, %f)\n", nearestPoint.x, nearestPoint.y, nearestPoint.z);


    // *********** free *************
    free_HashNode(adjacent_regions);
    free_HashNode(boundary_points);
    free_HashNode(voronoi_cells);
    free_GroupPoint(m_groupPoint.next);
    free_GroupPoint(seed_points);
}

int main()
{
    // test_0();
    test_1();

    return 0;
}