#include <nan.h>
#include <functional>
#include <iostream>

#include "exchange.h"
#include "fire.h"

using namespace Nan;
using namespace std;
using namespace v8;

class GPWorker : public AsyncWorker {
    public:
        GPWorker(Callback *callback, int under)
        : AsyncWorker(callback), under(under), primes(0) {}

        ~GPWorker() {}


        void Execute () {
            Exchange x(
                [&](void * data) {
                    primes.push_back(*((int *) data));
                }
            );

            startGP(under, (void*)&x);
        }

        // We have the results, and we're back in the event loop.
        void HandleOKCallback () {
            Nan:: HandleScope scope;

            v8::Local<v8::Array> results = New<v8::Array>(primes.size());
            int i = 0;
            for_each(primes.begin(), primes.end(),
                [&](int value) {
                    Nan::Set(results, i, New<v8::Number>(value));
                    i++;
            });


            Local<Value> argv[] = {
                Null(),
                results
            };

            callback->Call(2, argv);

        }
    private:
        int under;
        vector<int> primes;
};

// Asynchronous access to the `getPrimes()` function
NAN_METHOD(DoGP) {
    int under = To<int>(info[0]).FromJust();
    Callback *callback = new Callback(info[1].As<Function>());

    AsyncQueueWorker(new GPWorker(callback, under));
}


NAN_MODULE_INIT(Init) {
    Nan::Set(target, New<String>("getPrimes").ToLocalChecked(),
        GetFunction(New<FunctionTemplate>(DoGP)).ToLocalChecked());
}

NODE_MODULE(addon, Init)
