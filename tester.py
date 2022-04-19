import json
from os import path
import subprocess
import tempfile
import argparse
import time
from enum import Enum


def info(msg):
    print(f'[INFO] {msg}')


def warn(msg):
    print(f'\033[93m[WARNING] {msg}\033[0m')


def fail(msg):
    print(f'\033[91m[FAILURE] {msg}\033[0m')


def success(msg):
    print(f'\033[92m[SUCCESS] {msg}\033[0m')


def report_failure(status, stdout, stderr, expected_stdout, should_fail, goal):
    fail(f'Test failed with exit status {status}')
    if should_fail:
        info('Note that this test expects the program to fail and exit with a non-zero status code')
    info(f'Goal: {goal}')
    info(f'Stdout:\n{stdout}')
    info(f'Expected stdout:\n{expected_stdout}')
    if stderr != '':
        info(f'Stderr:\n{stderr}')


class TestType(Enum):
    C = 1
    Python = 2


def run_test(test_type, idx, file_contents, expected_output, should_fail, goal, k=None):
    info(f'Running test {idx}...')
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.csv') as tmp:
        tmp.write(file_contents)
        tmp.flush()

        start_time = time.time()
        if test_type == TestType.C:
            assert k is None
            result = subprocess.run(["./spkmeans", goal, tmp.name],
                                    capture_output=True, timeout=30, text=True, encoding='ascii')
        elif test_type == TestType.Python:
            assert k is not None
            result = subprocess.run(["python3", "spkmeans.py", str(k), goal, tmp.name],
                                    capture_output=True, timeout=30, text=True, encoding='ascii')
        else:
            print(f"Unknown test type {test_type}")
            exit(1)
        total_time = time.time() - start_time

        if should_fail:
            if result.returncode == 0 or result.stdout != expected_output:
                report_failure(result.returncode, result.stdout, result.stderr, expected_output, should_fail, goal)
                return False
            else:
                success(f'Test {idx} completed successfully in {total_time:.3f} seconds')
                return True
        else:
            if result.returncode != 0 or result.stdout != expected_output or result.stderr != '':
                report_failure(result.returncode, result.stdout, result.stderr, expected_output, should_fail, goal)
                return False
            else:
                success(f'Test {idx} completed successfully in {total_time:.3f} seconds')
                return True


def record_test(test_type, idx, file_contents, should_fail, goal, k=None):
    info(f'Running test {idx}...')
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.csv') as tmp:
        tmp.write(file_contents)
        tmp.flush()

        if test_type == TestType.C:
            assert k is None
            result = subprocess.run(["./spkmeans", goal, tmp.name],
                                    capture_output=True, timeout=30, text=True, encoding='ascii')
        elif test_type == TestType.Python:
            assert k is not None
            result = subprocess.run(["python3", "spkmeans.py", str(k), goal, tmp.name],
                                    capture_output=True, timeout=30, text=True, encoding='ascii')
        else:
            print(f"Unknown test type {test_type}")
            exit(1)

        if should_fail:
            if result.returncode == 0:
                fail('Zero return code found for should-fail test, skipping recording of this test')
            else:
                success(f'Test {idx} recorded successfully')
                return result.stdout
        else:
            if result.returncode != 0 or result.stderr != '':
                fail('Non-zero return code or errors found, skipping recording of this test')
            else:
                success(f'Test {idx} recorded successfully')
                return result.stdout


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--record', action='store_true', help='Record the current outputs as the test\'s expected output')
    args = parser.parse_args()

    with open('tests.json', 'r') as f:
        tests = json.loads(f.read())

    c_tests = tests['c']
    py_tests = tests['python']

    run_c_tests = path.isfile('spkmeans')
    run_py_tests = path.isfile('spkmeans.py')

    all_tests_passed = True
    all_recordings_made = True

    if run_c_tests:
        if not args.record:
            info('Running C tests:')
            for i, test in enumerate(c_tests):
                all_tests_passed = run_test(TestType.C, i, test['input'], test['expected_stdout'], test['should_fail'], test['goal']) and all_tests_passed
        else:
            info('Recording C tests:')
            for i, test in enumerate(c_tests):
                res = record_test(TestType.C, i, test['input'], test['should_fail'], test['goal'])
                if res is not None:
                    test['expected_stdout'] = res
                else:
                    all_recordings_made = False
    else:
        warn('C executable (./spkmeans) missing, skipping C tests...')

    if run_py_tests:
        if not args.record:
            info('Running Python tests:')
            for i, test in enumerate(py_tests):
                all_tests_passed = run_test(TestType.Python, i, test['input'], test['expected_stdout'], test['should_fail'], test['goal'], test['k']) and all_tests_passed
        else:
            info('Recording Python tests:')
            for i, test in enumerate(py_tests):
                res = record_test(TestType.Python, i, test['input'], test['should_fail'], test['goal'], test['k'])
                if res is not None:
                    test['expected_stdout'] = res
                else:
                    all_recordings_made = False
    else:
        warn('Python script (spkmeans.py) missing, skipping Python tests...')


    if args.record:
        if not all_recordings_made:
            fail('Not all recordings were successful! Scroll up for details')
        with open('tests.json', 'w') as f:
            f.write(json.dumps(tests, indent=4, sort_keys=True))
    else:
        if not all_tests_passed:
            fail('Not all tests were successful! Scroll up for details')



if __name__ == '__main__':
    main()