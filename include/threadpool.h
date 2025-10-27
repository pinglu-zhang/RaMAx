#pragma once
#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>
#include <type_traits>

// modified from https://github.com/progschj/ThreadPool.git
class ThreadPool {
public:
    ThreadPool(size_t);
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args)
        -> std::future<typename std::invoke_result<F, Args...>::type>;
    void waitAllTasksDone();
    ~ThreadPool();
private:
    // worker threads
    std::vector<std::thread> workers;
    // task queue
    std::queue<std::function<void()>> tasks;

    // synchronization
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;

    // bookkeeping for "all tasks done"
    std::condition_variable all_tasks_done;
    int tasks_count = 0;
};

// the constructor just launches some amount of workers
inline ThreadPool::ThreadPool(size_t threads)
    : stop(false)
{
    for (size_t i = 0; i < threads; ++i)
        workers.emplace_back(
            [this]
            {
                for (;;)
                {
                    std::function<void()> task;

                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock,
                            [this] { return this->stop || !this->tasks.empty(); });

                        // 如果线程池被停止并且任务队列也空了，就退出线程循环
                        if (this->stop && this->tasks.empty())
                            return;

                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }

                    // 关键修复：
                    // 防止用户任务抛异常导致整个 worker 线程直接死亡，
                    // 让线程池能长期保持满线程工作能力。
                    try {
                        task();
                    } catch (const std::exception &e) {
                        // 这里可以选择记录日志，例如：
                        // spdlog::error("ThreadPool worker caught std::exception: {}", e.what());
                        // 为了保持 header 干净，不引入 spdlog。
                    } catch (...) {
                        // spdlog::error("ThreadPool worker caught unknown exception");
                    }
                }
            }
        );
}

// add new work item to the pool
template<class F, class... Args>
auto ThreadPool::enqueue(F&& f, Args&&... args)
    -> std::future<typename std::invoke_result<F, Args...>::type>
{
    using return_type = typename std::invoke_result<F, Args...>::type;

    // 将任务绑定并封装成 packaged_task，这样我们可以拿到 future
    auto task = std::make_shared<std::packaged_task<return_type()>>(
        std::bind(std::forward<F>(f), std::forward<Args>(args)...)
    );

    std::future<return_type> res = task->get_future();
    {
        std::unique_lock<std::mutex> lock(queue_mutex);

        // 不允许在 stop 后继续 enqueue
        if (stop)
            throw std::runtime_error("enqueue on stopped ThreadPool");

        tasks_count++;

        // 这里放进队列的实际是一个无参 void() 包装
        tasks.emplace([task, this]() {
            // 运行用户任务
            (*task)();

            // 任务完成后更新计数并可能唤醒等待者
            std::unique_lock<std::mutex> lock(this->queue_mutex);
            tasks_count--;
            if (tasks_count == 0 && this->tasks.empty()) {
                all_tasks_done.notify_one();
            }
        });
    }
    condition.notify_one();
    return res;
}

inline void ThreadPool::waitAllTasksDone() {
    std::unique_lock<std::mutex> lock(queue_mutex);
    all_tasks_done.wait(lock, [this] {
        return tasks.empty() && tasks_count == 0;
    });
}

// the destructor joins all threads
inline ThreadPool::~ThreadPool()
{
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }

    // 通知所有 worker 线程退出等待
    condition.notify_all();

    // 回收所有 worker
    for (std::thread &worker : workers)
        worker.join();
}
