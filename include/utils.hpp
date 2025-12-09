#pragma once

#include <csignal>
#include <atomic>

std::atomic<bool> interrupt_requested(false);
void signal_handler(int signal) {
    if (signal == SIGINT) {
              interrupt_requested = true;
    }
}