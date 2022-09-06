#pragma once

#include <ve_offload.h>
// realpath関数を使用するため定数PATH_MAXを得る
#include <limits.h>
// dirname関数を使用するため libgen.hをインクルード
#include <libgen.h>
#include "time.h"

#if defined(_SXDEBUG)
#define VEO_DEBUG_PRINT(a,b) { std::cout<<" function: "<<__FUNCTION__<<" "<<(a)<<"="<<(b)<<std::endl; }
#define VEO_DEBUG_PRINT2(a,b,c,d) { std::cout<<" function: "<<__FUNCTION__<<" "<<(a)<<"="<<(b)<<" "<<(c)<<"="<<(d)<<std::endl; }
#else
#define VEO_DEBUG_PRINT(a,b)
#define VEO_DEBUG_PRINT2(a,b,c,d)
#endif

#if !defined(_SXDEBUG)
#define VEO_NO_PROFILE
#endif

#if !defined(VEO_NO_PROFILE)
#define VEO_USE_PROFILE
#endif

#if defined(VEO_METHOD2)
#define VEO_PROC_EACH_VECORE
#endif

// VEO_PROC_EACH_VECOREが定義されていればVEの１コア毎にプロセスを作成する。
// そうでなければ1VE毎にプロセスを作成し、8つのコンテキストを生成する。
// 開発段階では前者を方式２，後者を方式３と呼んだ。
// デフォルトは方式３とする。
//方式３では非同期のバッファ書き込みの方が速い

#if defined(VEO_PROC_EACH_VECORE)
#define VE_BUFFER_SYNC_RW
#else
#define VE_BUFFER_ASYNC_RW
#endif

#if !defined(PREPARE_EPJ_IN_VH)
#define PREPARE_EPJ_IN_VE
#endif

#if !defined(PREPARE_SPJ_IN_VE)
#define PREPARE_SPJ_IN_VH
#endif

//
// 注意
// 環境変数 VE_OMP_NUM_THREADは 1 に設定する必要がある。
// FDPSではOpenMPのスレッドからユーザ定義相互作用関数を呼ぶので
// VE側はOpenMPのマルチスレッドにしない。
// 一方、VE_OMP_NUM_THREADが2以上だと、VEプロセス生成時に複数プロセスが生成されてしまう。
//
// --- Alternative VE Offloading  2.7.5 の Restrictionsより ---
// When veo_proc_create() is invoked, multiple threads for a OpenMP program are created on VE side in the default context. 
// If you do not use OpenMP, set the environment variable VE_OMP_NUM_THREADS=1.
//

namespace VEOffloading {

#if defined(VEO_USE_PROFILE)
	// 時間取得関数
	inline double GetVeoTime() {
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
                return MPI_Wtime();
#elif defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
                return omp_get_wtime();
#else
                return (double)clock() / CLOCKS_PER_SEC;
#endif //PARTICLE_SIMULATOR_MPI_PARALLEL
	}
#endif	// VEO_USE_PROFILE

#define N_VEO_MAX_VE   8
#define N_VE_NUM_CORES 8
#define N_VEO_MAX_THREAD (N_VEO_MAX_VE * N_VE_NUM_CORES)

	// 環境変数の取得で使用する。
	class EnvVal {
        public:
		// 環境変数から整数値を得る
                static int getInt(const char* pEnvValName, int iMinVal, int iMaxVal) {
                        const char *pValStr = getenv(pEnvValName);
                        int retVal = iMinVal;
                        if (pValStr != NULL) {
                                try {
                                        retVal = std::stoi(pValStr);
                                }
                                catch (const std::invalid_argument& e) {
                                }
                                catch (const std::out_of_range& e) {
                                }
                        }

                        if (retVal < iMinVal) retVal = iMinVal;
                        if (retVal > iMaxVal) retVal = iMaxVal;

                        return retVal;
                }

		// OMP_NUM_THREADSを得る
		static int getNumThreads() {
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
                        return getInt("OMP_NUM_THREADS", 1, N_VEO_MAX_THREAD);
#else
			return 1;
#endif
		}

		// VE_NODE_NUMBER=0-3 なら 0,1,2,3のような配列に展開して返す
                static std::vector<int> getIntRange(const char* pEnvValName, int iMinVal, int iMaxVal) {
                        std::vector<int> arr;

                        const char *pValStr = getenv(pEnvValName);
                        std::string valstr(pValStr != NULL ? pValStr : "");


                        int iMin = 0;
                        int iMax = 0;

                        if (valstr.length() == 3 && valstr[1] == '-') {
                                try {
                                        std::string sMin;
                                        sMin.push_back(valstr[0]);
                                        iMin = std::stoi(sMin);
                                }
                                catch (const std::invalid_argument& e) {
                                }
                                catch (const std::out_of_range& e) {
                                }
                                try {
                                        std::string sMax;
                                        sMax.push_back(valstr[2]);
                                        iMax = std::stoi(sMax);
                                }
                                catch (const std::invalid_argument& e) {
                                }
                                catch (const std::out_of_range& e) {
                                }

                                if (iMin > iMax ) {
                                        std::swap(iMin,iMax);
                                }


                                if (iMin < iMinVal) iMin = iMinVal;
                                if (iMax > iMaxVal) iMax = iMaxVal;

                                for(int i=iMin; i <= iMax; i++) arr.push_back(i);
                        }
                        else if(valstr.length() == 1) {
                                try {
                                        std::string sMin;
                                        sMin.push_back(valstr[0]);
                                        iMax = iMin = std::stoi(sMin);
                                }
                                catch (const std::invalid_argument& e) {
                                }
                                catch (const std::out_of_range& e) {
                                }

                                if (iMin < iMinVal) iMax = iMin = iMinVal;
                                if (iMax > iMaxVal) iMax = iMin = iMaxVal;

                                for(int i=iMin; i <= iMax; i++) arr.push_back(i);
                        }
                        else if(valstr[1] == ',') {
                                for (int i=0; i<valstr.length(); i += 2) {
                                        int iVal = -1;
                                        try {
                                         std::string sVal;
                                         sVal.push_back(valstr[i]);
                                         iVal = std::stoi(sVal);
                                        }
                                        catch (const std::invalid_argument& e) {
                                        }
                                        catch (const std::out_of_range& e) {
                                        }

                                        if (iVal>=iMinVal && iVal<=iMaxVal && find(arr.begin(),arr.end(),iVal) == arr.end()) {
                                                arr.push_back(iVal);
                                        }
                                }
                        }

                        //for(int i=0; i<arr.size(); i++) {
                        //       std::cout << pEnvValName << "[" << i << "]" << "=" << arr[i]<< std::endl;
                        //}

                        return arr;
                }
	};

	// VEバッファへのアクセス回数・アクセス時間を保持する
	class VeoBufferProfile {
	public:
		const char* buffer_name_;
		int64_t  write_count_ ;
		int64_t  read_count_;
		int64_t	 write_bytes_;  
		int64_t  read_bytes_;
		double   write_time_;
		double	 read_time_;

		int	 veNum_;
		uint64_t symNum_; 
		char	funcName_[256];

		VeoBufferProfile(const char *buffer_name) : buffer_name_(buffer_name), write_count_(0), read_count_(0), write_bytes_(0), read_bytes_(0), write_time_(0.0), read_time_(0.0) {}

		void output() const {
                        if (write_count_ > 0) {
                                fprintf(stdout,"rank=%d veNum_=%d.%d funcName_=%s %s write_count_=%ld write_bytes_=%ld write_time_=%lf ave_write_bytes=%ld\n",
                                        ParticleSimulator::Comm::getRank(), veNum_, symNum_, funcName_, buffer_name_,
                                        write_count_, write_bytes_, write_time_, write_bytes_/write_count_);
                        }
                        if (read_count_ > 0) {
                                fprintf(stdout,"rank=%d veNum_=%d.%d funcName_=%s %s read_count_=%ld read_bytes_=%ld read_time_=%lf ave_read_bytes=%ld\n",
                                        ParticleSimulator::Comm::getRank(), veNum_, symNum_, funcName_, buffer_name_,
                                        read_count_, read_bytes_, read_time_, read_bytes_/read_count_);
                        }
		}
	};

	// VEOffloading呼び出しの回数と実行時間を保持する
	class VeoCallProfile {
	private:
		int64_t call_count_;
		double  total_call_time_;
		double  total_core_call_time_;
		int64_t total_epi_;
		int64_t total_epj_;
	public:
		int	veNum_;
		uint64_t symNum_; 
		char	funcName_[256];

		VeoCallProfile() : symNum_(999), call_count_(0), total_call_time_(0.0), total_core_call_time_(0.0), total_epi_(0), total_epj_(0)  {}

		void addCallData(double call_time, double core_call_time, int32_t n_epi, int32_t n_epj) {
			call_count_ += 1;
			total_call_time_ += call_time;
			total_core_call_time_ += core_call_time;
			total_epi_ += n_epi;
			total_epj_ += n_epj;
		}

		void output() const {
                        if (call_count_ > 0) {
                                fprintf(stdout,"rank=%d veNum_=%d.%d funcName_=%s call_count_=%ld total_call_time_=%lf total_core_call_time_=%lf ave_epi_=%ld ave_epj_=%ld\n",
                                        ParticleSimulator::Comm::getRank(), veNum_, symNum_, funcName_,
                                        call_count_, total_call_time_, total_core_call_time_, total_epi_/call_count_, total_epj_/call_count_);
                        }
		}
	};

	class VeoArgs {
	private:
		struct veo_args *argp_;
		int arg_count_;
	public:
		// コンストラクタ
		VeoArgs() {
			argp_ = veo_args_alloc();
			arg_count_ = 0;
		}
		// デストラクタ
		~VeoArgs() {
			veo_args_free(argp_);
		}

		struct veo_args * getArgp() const { return argp_; }

		void clear() {
			veo_args_clear(argp_);
			arg_count_ = 0;
		}

		void add(int arg_val) {
			veo_args_set_i32(argp_, arg_count_, arg_val);
			arg_count_++;
		}

		void add(uint64_t arg_val) {
			veo_args_set_u64(argp_, arg_count_, arg_val);
			arg_count_++;
		}
		void add(double arg_val) {
			veo_args_set_double(argp_, arg_count_, arg_val);
			arg_count_++;
		}
	};  // VeoArgs

#if defined(VE_BUFFER_SYNC_RW)
        class VeBuffer {
        private:
                struct veo_proc_handle *proc_;
                uint64_t bufptr_;
                size_t   bufsize_;

#if defined(VEO_USE_PROFILE)
                VeoBufferProfile profile_;
#endif
        public:
#if defined(VEO_USE_PROFILE)
		// コンストラクタ
		VeBuffer(struct veo_proc_handle *proc, int veNum, const char* buffer_name, int symNum, const char*  psym_name) : proc_(proc), bufptr_(0), bufsize_(0), profile_(buffer_name) {
                        profile_.veNum_ = veNum;
                        profile_.symNum_ = symNum;
                        strcpy(profile_.funcName_, psym_name);
                }
#else
		// コンストラクタ
                VeBuffer(struct veo_proc_handle *proc) : proc_(proc), bufptr_(0), bufsize_(0) {}
#endif
		// デストラクタ
                ~VeBuffer() { freeVeBuffer(); }

                void allocVeBuffer(size_t allocsize) {
                        if (bufsize_ > 0 && bufsize_ < allocsize) {
                                freeVeBuffer();
                        }

                        if (bufsize_ == 0) {
                                int ret = veo_alloc_mem(proc_, &bufptr_, allocsize);
                                bufsize_ = allocsize;
                        }
                }

                void freeVeBuffer() {
                        if (bufsize_ > 0) {
                                int ret = veo_free_mem(proc_, bufptr_);
                                bufsize_ = 0;
                        }
                }

                uint64_t getBufptr() const { return bufptr_; }

                int write(const void *buffer, size_t wsize) {
			// 今より大きなVEメモリが必要な場合には一旦割り当て解除
                        if (bufsize_ > 0 && bufsize_ < wsize) {
                                freeVeBuffer();
                        }

			// VEメモリが未割当なら割り当てを実施
                        if (bufsize_ == 0) {
                                allocVeBuffer(wsize);
                        }

#if defined(VEO_USE_PROFILE)
                        profile_.write_count_ += 1;
                        profile_.write_bytes_ += wsize;
                        double offset_time = GetVeoTime();
#endif
			// VEメモリに書き込み実行
                        int ret = veo_write_mem(proc_, bufptr_, buffer, wsize);
#if defined(VEO_USE_PROFILE)
                        profile_.write_time_ += GetVeoTime() - offset_time;
#endif
                        return ret;
                }

                int read(void *buffer, size_t rsize) {
                        if (rsize > bufsize_) {
                                return -1;
                        }
                        if (buffer == NULL) {
                                return -1;
                        }

#if defined(VEO_USE_PROFILE)
                        profile_.read_count_ += 1;
                        profile_.read_bytes_ += rsize;
                        double offset_time = GetVeoTime();
#endif
			// VEメモリから読み込みを実行
                        int ret = veo_read_mem(proc_, buffer, bufptr_, rsize);
#if defined(VEO_USE_PROFILE)
                        profile_.read_time_ += GetVeoTime() - offset_time;
#endif
                        return ret;
                }

#if defined(VEO_USE_PROFILE)
                void outputProfile() { profile_.output(); }
#endif
	};

#else   // not VE_BUFFER_SYNC_RW
	class VeBuffer {
	private:
		struct veo_proc_handle *proc_;
		veo_thr_ctxt	*ctx_;
		uint64_t bufptr_;			
		size_t	 bufsize_;

#if defined(VEO_USE_PROFILE)
		VeoBufferProfile profile_;
#endif
	public:
#if defined(VEO_USE_PROFILE)
		// コンストラクタ
		VeBuffer(struct veo_proc_handle *proc, veo_thr_ctxt *ctx, int veNum, const char* buffer_name, int symNum, const char*  psym_name) : proc_(proc), ctx_(ctx), bufptr_(0), bufsize_(0), profile_(buffer_name) {
			profile_.veNum_ = veNum;
			profile_.symNum_ = symNum;
			strcpy(profile_.funcName_, psym_name);
		}
#else
		// コンストラクタ
		VeBuffer(struct veo_proc_handle *proc, veo_thr_ctxt *ctx) : proc_(proc), ctx_(ctx), bufptr_(0), bufsize_(0) {}
#endif
		// デストラクタ
		~VeBuffer() { freeVeBuffer(); }	

		void allocVeBuffer(size_t allocsize) {
			if (bufsize_ > 0 && bufsize_ < allocsize) {
				freeVeBuffer();
			}

			if (bufsize_ == 0) {
				int ret = veo_alloc_mem(proc_, &bufptr_, allocsize);
				bufsize_ = allocsize;
			}
		}

		void freeVeBuffer() {
			if (bufsize_ > 0) {
				int ret = veo_free_mem(proc_, bufptr_);
				bufsize_ = 0;
			}
		}

		uint64_t getBufptr() const { return bufptr_; }

		int write(const void *buffer, size_t wsize) {
			// 今より大きなVEメモリが必要な場合には一旦割り当て解除
			if (bufsize_ > 0 && bufsize_ < wsize) {
				freeVeBuffer();
			}

			// VEメモリが未割当なら割り当てを実施
			if (bufsize_ == 0) {
				allocVeBuffer(wsize);
			}

#if defined(VEO_USE_PROFILE)
			profile_.write_count_ += 1;
			profile_.write_bytes_ += wsize;
			double offset_time = GetVeoTime();
#endif
			// VEメモリに書き込み実行
                        uint64_t id = veo_async_write_mem(ctx_, bufptr_, buffer, wsize);
                        uint64_t retval;
                        int ret = veo_call_wait_result(ctx_, id, &retval);
#if defined(VEO_USE_PROFILE)
			profile_.write_time_ += GetVeoTime() - offset_time;
#endif
			return ret;
		} 	

		int read(void *buffer, size_t rsize) {
			if (rsize > bufsize_) {
				return -1;
			}
			if (buffer == NULL) {
				return -1;
			}

#if defined(VEO_USE_PROFILE)
			profile_.read_count_ += 1;
			profile_.read_bytes_ += rsize;
			double offset_time = GetVeoTime();
#endif
			// VEメモリから読み込みを実行
                        uint64_t id = veo_async_read_mem(ctx_, buffer, bufptr_, rsize);
                        uint64_t retval;
                        int ret = veo_call_wait_result(ctx_, id, &retval);
#if defined(VEO_USE_PROFILE)
			profile_.read_time_ += GetVeoTime() - offset_time;
#endif
			return ret;
		}

#if defined(VEO_USE_PROFILE)
		void outputProfile() { profile_.output(); }
#endif
	};
#endif // end  VE_BUFFER_SYNC_RW or not

	class VeoContext {
	private:
		struct veo_thr_ctxt *ctx_;

	public:
		// コンストラクタ
		VeoContext(struct veo_proc_handle *proc) {
        		ctx_ = veo_context_open(proc);
			VEO_DEBUG_PRINT2("veo_context_open() ctx_", ctx_, "veo_num_contexts()", veo_num_contexts(proc));
		}

		// デストラクタ
		~VeoContext() {
			VEO_DEBUG_PRINT("veo_context_close() ctx_", ctx_);
			int close_status = veo_context_close(ctx_);
		}

		struct veo_thr_ctxt *getContext() const { return ctx_; }
	};

	class VeoProcess {
		// VEOプロセスハンドル
		struct veo_proc_handle *proc_;
		// ユーザ定義ダイナミックリンクライブラリのハンドル
		uint64_t libHandle_;

		// VEの最大コア数分の配列を用意
		VeoContext* pVeoContext_[N_VE_NUM_CORES];

		// VE番号
		int veNum_;
		// コンテキスト数
		int nContexts_;

	public:
		// コンストラクタ
		VeoProcess(const char *libpath, int venode, int nContexts) : veNum_(venode), nContexts_(nContexts) {
			for (int i=0; i<N_VE_NUM_CORES; i++) {
				pVeoContext_[i] = NULL;
			}

 			proc_ = veo_proc_create(venode);
        		if (proc_ == NULL) {
                		perror("veo_proc_create");
                		exit(1);
        		}

			VEO_DEBUG_PRINT("veo_proc_create proc_",proc_);
			VEO_DEBUG_PRINT2("libpath",libpath,"nContexts_",nContexts_);

			ParticleSimulator::F64 time_offset = ParticleSimulator::GetWtime();
        		libHandle_ = veo_load_library(proc_, libpath);
			VEO_DEBUG_PRINT("libHandle_", libHandle_);
			VEO_DEBUG_PRINT2("veo_load_library libHandle_", libHandle_, "time", ParticleSimulator::GetWtime() - time_offset);

			VEO_DEBUG_PRINT2("VeoProcess proc_", proc_, "libHandle_", libHandle_);
			if (libHandle_ == 0) {
				exit(-1);
			}

#if defined(VEO_PROC_EACH_VECORE)
			pVeoContext_[0] = new VeoContext(proc_);
#else
			for (int i=0; i<nContexts_; i++ ) {
				pVeoContext_[i] = new VeoContext(proc_);
			}
#endif
		}

		// デストラクタ
		~VeoProcess() {
			for (int i=0; i<N_VE_NUM_CORES; i++) {
				if (pVeoContext_[i] != NULL) delete pVeoContext_[i];
			}

			
			VEO_DEBUG_PRINT2("veo_unload_library proc_", proc_, "libHandle_", libHandle_);
			veo_unload_library(proc_, libHandle_);
			VEO_DEBUG_PRINT("veo_proc_destroy proc_", proc_);
			veo_proc_destroy(proc_);
		}

		struct veo_proc_handle* getProcHandle() const { return proc_; }

		uint64_t getLibHandle() const { return libHandle_; }

		struct veo_thr_ctxt* getContext(int ith) const {
			return pVeoContext_[ith]->getContext();
		}

		uint64_t getSymbol(const char* name) {
			return veo_get_sym(proc_, libHandle_, name);
		}

		int getVeNumber() const { return veNum_; }

		int getNumberOfContexts() const { return nContexts_; }
	};

	class InnerVeoFuncEpEp {
	private:
		VeoProcess* pVeoProcess_;
		uint64_t prepare_symid_;
		uint64_t symid_;
		bool  bSuperParticleType_;
		VeBuffer* pSortedEpiBuf_;
		VeBuffer* pSortedEpjBuf_;
		VeoArgs* pVeoArgs_[N_VE_NUM_CORES];
		VeBuffer* pAdrEpjForForceBuf_[N_VE_NUM_CORES];

#if defined(VEO_USE_PROFILE)
                VeoCallProfile prepare_profile_[N_VE_NUM_CORES];
                VeoCallProfile profile_[N_VE_NUM_CORES];
#endif
	public:
#if defined(VEO_USE_PROFILE)
		// コンストラクタ
		InnerVeoFuncEpEp(VeoProcess* pVeoProcess, int symNum, const char* psym_prepare_name, const char*  psym_name, InnerVeoFuncEpEp *pEpEpFunc, VeoCallProfile& prepare_profile, VeoCallProfile& profile) : pVeoProcess_(pVeoProcess) {
#else
		// コンストラクタ
		InnerVeoFuncEpEp(VeoProcess* pVeoProcess, int symNum, const char* psym_prepare_name, const char*  psym_name, InnerVeoFuncEpEp *pEpEpFunc) : pVeoProcess_(pVeoProcess) {
#endif
			for (int i=0; i < N_VE_NUM_CORES ; i++) {
				pVeoArgs_[i] = NULL;
				pAdrEpjForForceBuf_[i] = NULL;
			}

			if (strcmp(psym_prepare_name,"dummy") != 0) {
				prepare_symid_ = pVeoProcess_->getSymbol(psym_prepare_name);
				VEO_DEBUG_PRINT2("psym_prepare_name", psym_prepare_name, "psym_name", psym_name);
			}
			symid_ = pVeoProcess_->getSymbol(psym_name);
			VEO_DEBUG_PRINT2("prepare_symid_", prepare_symid_, "symid_", symid_);
                        bSuperParticleType_ = (pEpEpFunc != NULL);
                        if (bSuperParticleType_) {
                                //pSortedForceBuf_ = pEpEpFunc->pSortedForceBuf_;
				//for (int i=0; i < N_VE_NUM_CORES ; i++) {
				//	pForceBuf_[i] = pEpEpFunc->pForceBuf_[i];
				//}
                                pSortedEpiBuf_ = pEpEpFunc->pSortedEpiBuf_;
                        }

#if defined(VE_BUFFER_SYNC_RW)
#if defined(VEO_USE_PROFILE)
#if defined(VEO_PROC_EACH_VECORE)
			if (!bSuperParticleType_) {
				pSortedEpiBuf_ = new VeBuffer(pVeoProcess_->getProcHandle(),
						pVeoProcess->getVeNumber(), "sorted_epi", symNum, psym_name);
			}
			pSortedEpjBuf_ = new VeBuffer(pVeoProcess_->getProcHandle(), 
						pVeoProcess->getVeNumber(), "sorted_epj", symNum, psym_name);
#else
			if (!bSuperParticleType_) {
				pSortedEpiBuf_ = new VeBuffer(pVeoProcess_->getProcHandle(),
						pVeoProcess->getVeNumber(), "sorted_epi", 0, psym_name);
			}
			pSortedEpjBuf_ = new VeBuffer(pVeoProcess_->getProcHandle(), 
						pVeoProcess->getVeNumber(), "sorted_epj", 0, psym_name);
#endif
			for (int i=0; i < pVeoProcess_->getNumberOfContexts(); i++) {
				pVeoArgs_[i] = new VeoArgs();
#if defined(VEO_PROC_EACH_VECORE)
				pAdrEpjForForceBuf_[i] = new VeBuffer(pVeoProcess_->getProcHandle(), 
								pVeoProcess->getVeNumber(), "adr_epj_for_force", symNum, psym_name);
#else
				pAdrEpjForForceBuf_[i] = new VeBuffer(pVeoProcess_->getProcHandle(), 
								pVeoProcess->getVeNumber(), "adr_epj_for_force", i % N_VE_NUM_CORES, psym_name);
#endif

                        	strcpy(prepare_profile_[i].funcName_, psym_prepare_name);
                        	strcpy(profile_[i].funcName_, psym_name);
                        	prepare_profile_[i].veNum_ = profile_[i].veNum_ = pVeoProcess->getVeNumber();
#if defined(VEO_PROC_EACH_VECORE)
                        	prepare_profile_[i].symNum_ = profile_[i].symNum_ = symNum;
#else
                        	prepare_profile_[i].symNum_ = profile_[i].symNum_ = i % N_VE_NUM_CORES;
#endif
			}
#else	// not VEO_USE_PROFILE
                        if (!bSuperParticleType_) {
                                pSortedEpiBuf_ = new VeBuffer(pVeoProcess_->getProcHandle());
                        }
			pSortedEpjBuf_ = new VeBuffer(pVeoProcess_->getProcHandle());
			for (int i=0; i < pVeoProcess_->getNumberOfContexts(); i++) {
				pVeoArgs_[i] = new VeoArgs();
				pAdrEpjForForceBuf_[i] = new VeBuffer(pVeoProcess_->getProcHandle());
			}
#endif
#else // not VE_BUFFER_SYNC_RW)
#if defined(VEO_USE_PROFILE)
			if (!bSuperParticleType_) {
				pSortedEpiBuf_ = new VeBuffer(pVeoProcess_->getProcHandle(), pVeoProcess_->getContext(0),
						pVeoProcess->getVeNumber(), "sorted_epi", 0, psym_name);
			}
			pSortedEpjBuf_ = new VeBuffer(pVeoProcess_->getProcHandle(), pVeoProcess_->getContext(0),
						pVeoProcess->getVeNumber(), "sorted_epj", 0, psym_name);
			for (int i=0; i < pVeoProcess_->getNumberOfContexts(); i++) {
				pVeoArgs_[i] = new VeoArgs();
				pAdrEpjForForceBuf_[i] = new VeBuffer(pVeoProcess_->getProcHandle(), pVeoProcess_->getContext(i),
								pVeoProcess->getVeNumber(), "adr_epj_for_force", i % N_VE_NUM_CORES, psym_name);

                        	strcpy(prepare_profile_[i].funcName_, psym_prepare_name);
                        	strcpy(profile_[i].funcName_, psym_name);
                        	prepare_profile_[i].veNum_ = profile_[i].veNum_ = pVeoProcess->getVeNumber();
                        	prepare_profile_[i].symNum_ = profile_[i].symNum_ = i % N_VE_NUM_CORES;
			}
#else
                        if (!bSuperParticleType_) {
                                pSortedEpiBuf_ = new VeBuffer(pVeoProcess_->getProcHandle(), pVeoProcess_->getContext(0));
                        }
			pSortedEpjBuf_ = new VeBuffer(pVeoProcess_->getProcHandle(), pVeoProcess_->getContext(0));
			for (int i=0; i < pVeoProcess_->getNumberOfContexts(); i++) {
				pVeoArgs_[i] = new VeoArgs();
				pAdrEpjForForceBuf_[i] = new VeBuffer(pVeoProcess_->getProcHandle(), pVeoProcess_->getContext(i));
			}
#endif
#endif  // end VE_BUFFER_SYNC_RW
		}

		// デストラクタ
		~InnerVeoFuncEpEp() {
#if defined(VEO_USE_PROFILE)
			// profileを出力する
			if (!bSuperParticleType_) {
				pSortedEpiBuf_->outputProfile();
			}
			pSortedEpjBuf_->outputProfile();
			for (int i=0; i < pVeoProcess_->getNumberOfContexts(); i++) {
				if (pAdrEpjForForceBuf_[i] != NULL) pAdrEpjForForceBuf_[i]->outputProfile();
                        	prepare_profile_[i].output();
                        	profile_[i].output();
			}
#endif
			if (!bSuperParticleType_) {
				delete pSortedEpiBuf_;
			}
			delete pSortedEpjBuf_;

			for (int i=0; i < N_VE_NUM_CORES ; i++) {
				if (pAdrEpjForForceBuf_[i] != NULL) delete pAdrEpjForForceBuf_[i];
				if (pVeoArgs_[i] != NULL) delete pVeoArgs_[i];
			}
		}

		VeBuffer* createVeBuffer(int iCtx) {
#if defined(VE_BUFFER_SYNC_RW)
#if defined(VEO_USE_PROFILE)
			return new VeBuffer(pVeoProcess_->getProcHandle(), 
								pVeoProcess_->getVeNumber(), "receive_force", iCtx, "func1");
#else
			return new VeBuffer(pVeoProcess_->getProcHandle());
#endif
#else // not VE_BUFFER_SYNC_RW
#if defined(VEO_USE_PROFILE)
			return new VeBuffer(pVeoProcess_->getProcHandle(), pVeoProcess_->getContext(iCtx),
								pVeoProcess_->getVeNumber(), "receive_force", iCtx, "func1");
#else
			return new VeBuffer(pVeoProcess_->getProcHandle(), pVeoProcess_->getContext(iCtx));
#endif
#endif
		}

		void downloadSortedEpEp(const void* pSortedEpiBuf, size_t epi_size, const void* pSortedEpjBuf, size_t epj_size) {
			if (!bSuperParticleType_) {
				pSortedEpiBuf_->write(pSortedEpiBuf, epi_size);
				//VEO_DEBUG_PRINT2("download pSortedEpiBuf_", pSortedEpiBuf_, "epi_size", epi_size);
			}
#if defined(PREPARE_EPJ_IN_VE)
			if (epj_size > 0) {
				pSortedEpjBuf_->write(pSortedEpjBuf, epj_size);
				//VEO_DEBUG_PRINT2("download pSortedEpjBuf_", pSortedEpjBuf_, "epj_size", epj_size);
			}
#endif
		}

		void downloadSortedSp(const void* pSortedEpjBuf, size_t epj_size) {
#if defined(PREPARE_SPJ_IN_VE)  
			if (epj_size > 0) {
				pSortedEpjBuf_->write(pSortedEpjBuf, epj_size);
				//VEO_DEBUG_PRINT2("download pSortedEpjBuf_", pSortedEpjBuf_, "spj_size", epj_size);
			}
#endif
		}

#if defined(VEO_PROC_EACH_VECORE)	// sync call
#if defined(PREPARE_EPJ_IN_VE)
		int callFuncEpEp(int ith, size_t offset_epi_size, int n_epi,
				const ParticleSimulator::S32* padr_epj_force_buf, size_t epj_size, int n_epj,
				VeBuffer *pVeEpjForForceBuf,
				VeBuffer *pVeForceBuf, void *pForceBuf, size_t force_size, const double *params = NULL,  int params_count = 0, bool bLongType = false) {
#if defined(VEO_USE_PROFILE)
                        double start_offset_time = GetVeoTime();
#endif
                        pAdrEpjForForceBuf_[0]->write(padr_epj_force_buf, sizeof(ParticleSimulator::S32) * n_epj);
                        pVeEpjForForceBuf->allocVeBuffer(epj_size);

                        pVeoArgs_[0]->clear();
                        pVeoArgs_[0]->add(pSortedEpjBuf_->getBufptr());
                        pVeoArgs_[0]->add(pAdrEpjForForceBuf_[0]->getBufptr());
                        pVeoArgs_[0]->add(n_epj);
                        pVeoArgs_[0]->add(pVeEpjForForceBuf->getBufptr());

#if defined(VEO_USE_PROFILE)
                        double prepare_offset_time = GetVeoTime();
#endif
                        uint64_t retval;
			int ret = veo_call_sync(pVeoProcess_->getProcHandle(), prepare_symid_, pVeoArgs_[0]->getArgp(), &retval);
#if defined(VEO_USE_PROFILE)
                        prepare_profile_[ith].addCallData(GetVeoTime() - start_offset_time, GetVeoTime() - prepare_offset_time, 0, 0) ;
                        double start_offset_time2 = GetVeoTime();
#endif
                        pVeForceBuf->allocVeBuffer(force_size);

                        pVeoArgs_[0]->clear();
                        pVeoArgs_[0]->add(pSortedEpiBuf_->getBufptr()+offset_epi_size);
                        pVeoArgs_[0]->add(n_epi);
                        pVeoArgs_[0]->add(pVeEpjForForceBuf->getBufptr());
                        pVeoArgs_[0]->add(n_epj);
                        pVeoArgs_[0]->add(pVeForceBuf->getBufptr());
                        if (params != NULL && params_count > 0) {
				for (int i=0; i<params_count; i++) {
                        		pVeoArgs_[0]->add(params[i]);
				}
			}

#if defined(VEO_USE_PROFILE)
                        double offset_time = GetVeoTime();
#endif
			ret = veo_call_sync(pVeoProcess_->getProcHandle(), symid_, pVeoArgs_[0]->getArgp(), &retval);
#if defined(VEO_USE_PROFILE)
                        profile_[ith].addCallData(GetVeoTime() - start_offset_time, GetVeoTime() - offset_time, n_epi, n_epj) ;
#endif

			if (!bLongType) {
                        	ret = pVeForceBuf->read(pForceBuf, force_size);
			}

                        return ret;
		}
#else // not defined PREPARE_EPJ_IN_VE   prepare epj in VH
                int callFuncEpEp(int ith, size_t offset_epi_size, int n_epi,
                                const void* epj_for_force_buf, size_t epj_size, int n_epj,
				VeBuffer *pVeEpjForForceBuf,
                                VeBuffer *pVeForceBuf, void *pForceBuf, size_t force_size, const double *params = NULL,  int params_count = 0, bool bLongType = false) {
#if defined(VEO_USE_PROFILE)
                        double start_offset_time = GetVeoTime();
#endif
                        pVeEpjForForceBuf->write(epj_for_force_buf, epj_size);
                        pVeForceBuf->allocVeBuffer(force_size);

                        pVeoArgs_[0]->clear();
                        pVeoArgs_[0]->add(pSortedEpiBuf_->getBufptr()+offset_epi_size);
                        pVeoArgs_[0]->add(n_epi);
                        pVeoArgs_[0]->add(pVeEpjForForceBuf->getBufptr());
                        pVeoArgs_[0]->add(n_epj);
                        pVeoArgs_[0]->add(pVeForceBuf->getBufptr());
                        if (params != NULL && params_count > 0) {
                                for (int i=0; i<params_count; i++) {
                                        pVeoArgs_[0]->add(params[i]);
                                }
                        }

#if defined(VEO_USE_PROFILE)
                        double offset_time = GetVeoTime();
#endif
                        uint64_t retval;
			int ret = veo_call_sync(pVeoProcess_->getProcHandle(), symid_, pVeoArgs_[0]->getArgp(), &retval);

#if defined(VEO_USE_PROFILE)
                        profile_[ith].addCallData(GetVeoTime() - start_offset_time, GetVeoTime() - offset_time, n_epi, n_epj) ;
#endif

			if (!bLongType) {
                        	ret = pVeForceBuf->read(pForceBuf, force_size);
			}

                        return ret;
                }
#endif	// end PREPARE_EPJ_IN_VE


#else // not VEO_PROC_EACH_VECORE	// async call

#if defined(PREPARE_EPJ_IN_VE)
		int callFuncEpEp(int ith, size_t offset_epi_size, int n_epi,
				const ParticleSimulator::S32* padr_epj_force_buf, size_t epj_size, int n_epj,
				VeBuffer *pVeEpjForForceBuf,
				VeBuffer *pVeForceBuf, void *pForceBuf, size_t force_size, const double *params = NULL,  int params_count = 0, bool bLongType = false) {
#if defined(VEO_USE_PROFILE)
                        double start_offset_time = GetVeoTime();
#endif
                        pAdrEpjForForceBuf_[ith]->write(padr_epj_force_buf, sizeof(ParticleSimulator::S32) * n_epj);
                        pVeEpjForForceBuf->allocVeBuffer(epj_size);

                        pVeoArgs_[ith]->clear();
                        pVeoArgs_[ith]->add(pSortedEpjBuf_->getBufptr());
                        pVeoArgs_[ith]->add(pAdrEpjForForceBuf_[ith]->getBufptr());
                        pVeoArgs_[ith]->add(n_epj);
                        pVeoArgs_[ith]->add(pVeEpjForForceBuf->getBufptr());

#if defined(VEO_USE_PROFILE)
                        double prepare_offset_time = GetVeoTime();
#endif
                        uint64_t retval;
                       	uint64_t id = veo_call_async(pVeoProcess_->getContext(ith), prepare_symid_, pVeoArgs_[ith]->getArgp());
                        veo_call_wait_result(pVeoProcess_->getContext(ith), id, &retval);

#if defined(VEO_USE_PROFILE)
                        prepare_profile_[ith].addCallData(GetVeoTime() - start_offset_time, GetVeoTime() - prepare_offset_time, 0, 0) ;
                        double start_offset_time2 = GetVeoTime();
#endif
                        pVeForceBuf->allocVeBuffer(force_size);

                        pVeoArgs_[ith]->clear();
                        pVeoArgs_[ith]->add(pSortedEpiBuf_->getBufptr()+offset_epi_size);
                        pVeoArgs_[ith]->add(n_epi);
                        pVeoArgs_[ith]->add(pVeEpjForForceBuf->getBufptr());
                        pVeoArgs_[ith]->add(n_epj);
                        pVeoArgs_[ith]->add(pVeForceBuf->getBufptr());
                        if (params != NULL && params_count > 0) {
				for (int i=0; i<params_count; i++) {
                        		pVeoArgs_[ith]->add(params[i]);
				}
			}

#if defined(VEO_USE_PROFILE)
                        double offset_time = GetVeoTime();
#endif
                       	id = veo_call_async(pVeoProcess_->getContext(ith), symid_, pVeoArgs_[ith]->getArgp());
                        veo_call_wait_result(pVeoProcess_->getContext(ith), id, &retval);

#if defined(VEO_USE_PROFILE)
                        profile_[ith].addCallData(GetVeoTime() - start_offset_time2, GetVeoTime() - offset_time, n_epi, n_epj) ;
#endif

                        int ret = 0;
			if (!bLongType) {
                        	ret = pVeForceBuf->read(pForceBuf, force_size);
			}

                        return ret;
		}
#else // not defined PREPARE_EPJ_IN_VE  prepare epj in VH
                int callFuncEpEp(int ith, size_t offset_epi_size, int n_epi,
                                const void* epj_for_force_buf, size_t epj_size, int n_epj,
				VeBuffer *pVeEpjForForceBuf,
                                VeBuffer *pVeForceBuf, void *pForceBuf, size_t force_size, const double *params = NULL,  int params_count = 0, bool bLongType = false) {
#if defined(VEO_USE_PROFILE)
                        double start_offset_time = GetVeoTime();
#endif
                        pVeEpjForForceBuf->write(epj_for_force_buf, epj_size);
                        pVeForceBuf->allocVeBuffer(force_size);

                        pVeoArgs_[ith]->clear();
                        pVeoArgs_[ith]->add(pSortedEpiBuf_->getBufptr()+offset_epi_size);
                        pVeoArgs_[ith]->add(n_epi);
                        pVeoArgs_[ith]->add(pVeEpjForForceBuf->getBufptr());
                        pVeoArgs_[ith]->add(n_epj);
                        pVeoArgs_[ith]->add(pVeForceBuf->getBufptr());
                        if (params != NULL && params_count > 0) {
                                for (int i=0; i<params_count; i++) {
                                        pVeoArgs_[ith]->add(params[i]);
                                }
                        }

#if defined(VEO_USE_PROFILE)
                        double offset_time = GetVeoTime();
#endif
                        uint64_t retval;
                        uint64_t id = veo_call_async(pVeoProcess_->getContext(ith), symid_, pVeoArgs_[ith]->getArgp());
                        veo_call_wait_result(pVeoProcess_->getContext(ith), id, &retval);

#if defined(VEO_USE_PROFILE)
                        profile_[ith].addCallData(GetVeoTime() - start_offset_time, GetVeoTime() - offset_time, n_epi, n_epj) ;
#endif

		
                        int ret = 0;
			if (!bLongType) {
                        	ret = pVeForceBuf->read(pForceBuf, force_size);
			}
                        return ret;
                }
#endif
#endif

#if defined(VEO_PROC_EACH_VECORE)	// sync call
#if defined(PREPARE_SPJ_IN_VE)
		int callFuncEpSp(int ith, size_t offset_epi_size, int n_epi,
				ParticleSimulator::S32* padr_spj_for_forcebuf, size_t spj_size, int n_spj,
				VeBuffer *pVeEpjForForceBuf,
				VeBuffer *pVeForceBuf, void *pForceBuf, size_t force_size, const double *params = NULL, int params_count = 0) {
#if defined(VEO_USE_PROFILE)
                        double start_offset_time = GetVeoTime();
#endif
                        pAdrEpjForForceBuf_[0]->write(padr_spj_for_forcebuf, sizeof(ParticleSimulator::S32) * n_spj);
                        pVeEpjForForceBuf->allocVeBuffer(spj_size);

                        pVeoArgs_[0]->clear();
                        pVeoArgs_[0]->add(pSortedEpjBuf_->getBufptr());
                        pVeoArgs_[0]->add(pAdrEpjForForceBuf_[0]->getBufptr());
                        pVeoArgs_[0]->add(n_spj);
                        pVeoArgs_[0]->add(pVeEpjForForceBuf->getBufptr());

#if defined(VEO_USE_PROFILE)
                        double prepare_offset_time = GetVeoTime();
#endif
                        uint64_t retval;
			int ret = veo_call_sync(pVeoProcess_->getProcHandle(), prepare_symid_, pVeoArgs_[0]->getArgp(), &retval);
#if defined(VEO_USE_PROFILE)
                        prepare_profile_[ith].addCallData(GetVeoTime() - start_offset_time, GetVeoTime() - prepare_offset_time, n_epi, n_spj) ;
                        double start_offset_time2 = GetVeoTime();
#endif

                        pVeoArgs_[0]->clear();
                        pVeoArgs_[0]->add(pSortedEpiBuf_->getBufptr()+offset_epi_size);
                        pVeoArgs_[0]->add(n_epi);
                        pVeoArgs_[0]->add(pVeEpjForForceBuf->getBufptr());
                        pVeoArgs_[0]->add(n_spj);
                        pVeoArgs_[0]->add(pVeForceBuf->getBufptr());
                        if (params != NULL && params_count > 0) {
				for (int i=0; i<params_count; i++) {
                        		pVeoArgs_[0]->add(params[i]);
				}
			}

#if defined(VEO_USE_PROFILE)
                        double offset_time = GetVeoTime();
#endif
			ret = veo_call_sync(pVeoProcess_->getProcHandle(), symid_, pVeoArgs_[0]->getArgp(), &retval);

#if defined(VEO_USE_PROFILE)
                        profile_[ith].addCallData(GetVeoTime() - start_offset_time2, GetVeoTime() - offset_time, n_epi, n_spj) ;
#endif
                        ret = pVeForceBuf->read(pForceBuf, force_size);
                        return ret;
		}
#else	// not PREPARE_SPJ_IN_VE
		int callFuncEpSp(int ith, size_t offset_epi_size, int n_epi,
				const void* spj_for_forcebuf, size_t spj_size, int n_spj,
				VeBuffer *pVeEpjForForceBuf,
				VeBuffer *pVeForceBuf, void *pForceBuf, size_t force_size, const double *params = NULL, int params_count = 0) {
#if defined(VEO_USE_PROFILE)
                        double start_offset_time = GetVeoTime();
#endif
                        pVeEpjForForceBuf->write(spj_for_forcebuf, spj_size);

                        pVeoArgs_[0]->clear();
                        pVeoArgs_[0]->add(pSortedEpiBuf_->getBufptr()+offset_epi_size);
                        pVeoArgs_[0]->add(n_epi);
                        pVeoArgs_[0]->add(pVeEpjForForceBuf->getBufptr());
                        pVeoArgs_[0]->add(n_spj);
                        pVeoArgs_[0]->add(pVeForceBuf->getBufptr());
                        if (params != NULL && params_count > 0) {
				for (int i=0; i<params_count; i++) {
                        		pVeoArgs_[0]->add(params[i]);
				}
			}

#if defined(VEO_USE_PROFILE)
                        double offset_time = GetVeoTime();
#endif
                        uint64_t retval;
			int ret = veo_call_sync(pVeoProcess_->getProcHandle(), symid_, pVeoArgs_[0]->getArgp(), &retval);

#if defined(VEO_USE_PROFILE)
                        profile_[ith].addCallData(GetVeoTime() - start_offset_time, GetVeoTime() - offset_time, n_epi, n_spj) ;
#endif
                        ret = pVeForceBuf->read(pForceBuf, force_size);
                        return ret;
		}
#endif	// end PREPARE_SPJ_IN_VE
#else // not VEO_PROC_EACH_VECORE  async call
#if defined(PREPARE_SPJ_IN_VE)
		int callFuncEpSp(int ith, size_t offset_epi_size, int n_epi,
				ParticleSimulator::S32* padr_spj_for_forcebuf, size_t spj_size, int n_spj,
				VeBuffer *pVeEpjForForceBuf,
				VeBuffer *pVeForceBuf, void *pForceBuf, size_t force_size, const double *params = NULL, int params_count = 0) {
#if defined(VEO_USE_PROFILE)
                        double start_offset_time = GetVeoTime();
#endif
                        pAdrEpjForForceBuf_[ith]->write(padr_spj_for_forcebuf, sizeof(ParticleSimulator::S32) * n_spj);
                        pVeEpjForForceBuf->allocVeBuffer(spj_size);

                        pVeoArgs_[ith]->clear();
                        pVeoArgs_[ith]->add(pSortedEpjBuf_->getBufptr());
                        pVeoArgs_[ith]->add(pAdrEpjForForceBuf_[ith]->getBufptr());
                        pVeoArgs_[ith]->add(n_spj);
                        pVeoArgs_[ith]->add(pVeEpjForForceBuf->getBufptr());

#if defined(VEO_USE_PROFILE)
                        double prepare_offset_time = GetVeoTime();
#endif
                        uint64_t retval;
                       	uint64_t id = veo_call_async(pVeoProcess_->getContext(ith), prepare_symid_, pVeoArgs_[ith]->getArgp());
                        veo_call_wait_result(pVeoProcess_->getContext(ith), id, &retval);

#if defined(VEO_USE_PROFILE)
                        prepare_profile_[ith].addCallData(GetVeoTime() - start_offset_time, GetVeoTime() - prepare_offset_time, n_epi, n_spj) ;
                        double start_offset_time2 = GetVeoTime();
#endif
//ccccc
                        pVeoArgs_[ith]->clear();
                        pVeoArgs_[ith]->add(pSortedEpiBuf_->getBufptr()+offset_epi_size);
                        pVeoArgs_[ith]->add(n_epi);
                        pVeoArgs_[ith]->add(pVeEpjForForceBuf->getBufptr());
                        pVeoArgs_[ith]->add(n_spj);
                        pVeoArgs_[ith]->add(pVeForceBuf->getBufptr());
                        if (params != NULL && params_count > 0) {
				for (int i=0; i<params_count; i++) {
                        		pVeoArgs_[ith]->add(params[i]);
				}
			}

#if defined(VEO_USE_PROFILE)
                        double offset_time = GetVeoTime();
#endif
                       	id = veo_call_async(pVeoProcess_->getContext(ith), symid_, pVeoArgs_[ith]->getArgp());
                        veo_call_wait_result(pVeoProcess_->getContext(ith), id, &retval);

#if defined(VEO_USE_PROFILE)
                        profile_[ith].addCallData(GetVeoTime() - start_offset_time2, GetVeoTime() - offset_time, n_epi, n_spj) ;
#endif
                        int ret = pVeForceBuf->read(pForceBuf, force_size);
                        return ret;
		}
//bbbbb
#else	// not PREPARE_SPJ_IN_VE
		int callFuncEpSp(int ith, size_t offset_epi_size, int n_epi,
				const void* spj_for_forcebuf, size_t spj_size, int n_spj,
				VeBuffer *pVeEpjForForceBuf,
				VeBuffer *pVeForceBuf, void *pForceBuf, size_t force_size, const double *params = NULL, int params_count = 0) {
#if defined(VEO_USE_PROFILE)
                        double start_offset_time = GetVeoTime();
#endif
                        pVeEpjForForceBuf->write(spj_for_forcebuf, spj_size);

                        pVeoArgs_[ith]->clear();
                        pVeoArgs_[ith]->add(pSortedEpiBuf_->getBufptr()+offset_epi_size);
                        pVeoArgs_[ith]->add(n_epi);
                        pVeoArgs_[ith]->add(pVeEpjForForceBuf->getBufptr());
                        pVeoArgs_[ith]->add(n_spj);
                        pVeoArgs_[ith]->add(pVeForceBuf->getBufptr());
                        if (params != NULL && params_count > 0) {
				for (int i=0; i<params_count; i++) {
                        		pVeoArgs_[ith]->add(params[i]);
				}
			}

#if defined(VEO_USE_PROFILE)
                        double offset_time = GetVeoTime();
#endif
                        uint64_t retval;
                       	uint64_t id = veo_call_async(pVeoProcess_->getContext(ith), symid_, pVeoArgs_[ith]->getArgp());
                        veo_call_wait_result(pVeoProcess_->getContext(ith), id, &retval);

#if defined(VEO_USE_PROFILE)
                        profile_[ith].addCallData(GetVeoTime() - start_offset_time, GetVeoTime() - offset_time, n_epi, n_spj) ;
#endif
                        int ret = pVeForceBuf->read(pForceBuf, force_size);
                        return ret;
		}
#endif	// end PREPARE_SPJ_IN_VE
#endif	// end  VEO_PROC_EACH_VECORE
	};

	class VeoFuncEpEp {
	private:
#if defined(VEO_PROC_EACH_VECORE)
		InnerVeoFuncEpEp* pInnerVeoFuncEpEp_[N_VEO_MAX_THREAD];	
#else
		InnerVeoFuncEpEp* pInnerVeoFuncEpEp_[N_VEO_MAX_VE];	
#endif
		VeBuffer*	pEpjForForceBuf_[N_VEO_MAX_THREAD];
		VeBuffer*	pForceBuf_[N_VEO_MAX_THREAD];
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
                //omp_lock_t      omp_lock_[N_VEO_MAX_THREAD];
                omp_lock_t      *p_omp_lock_;
#endif
		int             n_veo_thread_count;

		bool 		bSuperParticleType_;

		bool 		bLongType_;	

#if defined(VEO_USE_PROFILE)
		VeoCallProfile prepare_profile_;
		VeoCallProfile profile_;
#endif

	public:
		// コンストラクタ
		VeoFuncEpEp(VeoProcess **ppVeoProcess, int symNum, const char* psym_prepare_name, const char*  psym_name, VeoFuncEpEp *pEpEpFunc) {
			bLongType_ = false;
                        n_veo_thread_count = 0;

#if defined(VEO_PROC_EACH_VECORE)
                        for (int i=0; i<N_VEO_MAX_THREAD; i++) {
                                pInnerVeoFuncEpEp_[i] = NULL;

                                if (ppVeoProcess[i] != NULL) {
                                        InnerVeoFuncEpEp *pInnerEpiFunc = NULL;
                                        if (pEpEpFunc != NULL) {
                                                pInnerEpiFunc = pEpEpFunc->pInnerVeoFuncEpEp_[i];
                                        }
#if defined(VEO_USE_PROFILE)
                                        pInnerVeoFuncEpEp_[i] = new InnerVeoFuncEpEp(ppVeoProcess[i], i % N_VE_NUM_CORES, psym_prepare_name, psym_name, pInnerEpiFunc, prepare_profile_, profile_);
#else
                                        pInnerVeoFuncEpEp_[i] = new InnerVeoFuncEpEp(ppVeoProcess[i], i % N_VE_NUM_CORES, psym_prepare_name, psym_name, pInnerEpiFunc);
#endif
                                        n_veo_thread_count += 1;
       		                 	//VEO_DEBUG_PRINT2("symNum", symNum, "n_veo_thread_count", n_veo_thread_count);
                                }
                        }
#else // not VEO_PROC_EACH_VECORE
			for (int i=0; i<N_VEO_MAX_VE; i++) {
				pInnerVeoFuncEpEp_[i] = NULL;

				if (ppVeoProcess[i] != NULL) {
                                        InnerVeoFuncEpEp *pInnerEpiFunc = NULL;
                                        if (pEpEpFunc != NULL) {
                                                pInnerEpiFunc = pEpEpFunc->pInnerVeoFuncEpEp_[i];
                                        }
#if defined(VEO_USE_PROFILE)
					pInnerVeoFuncEpEp_[i] = new InnerVeoFuncEpEp(ppVeoProcess[i], i % N_VE_NUM_CORES, psym_prepare_name, psym_name, pInnerEpiFunc, prepare_profile_, profile_);
#else
					pInnerVeoFuncEpEp_[i] = new InnerVeoFuncEpEp(ppVeoProcess[i], i % N_VE_NUM_CORES, psym_prepare_name, psym_name, pInnerEpiFunc);
#endif

	                        	n_veo_thread_count += N_VE_NUM_CORES;
       		                 	//VEO_DEBUG_PRINT2("symNum", symNum, "n_veo_thread_count", n_veo_thread_count);
				}
			}
#endif

/*
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
                        for (int i=0; i<n_veo_thread_count; i++) {
                                omp_init_lock(&omp_lock_[i]);
			}
#endif
*/
			bSuperParticleType_ = false;
                        for (int i=0; i<N_VEO_MAX_THREAD; i++) {
				pEpjForForceBuf_[i] = NULL;
				pForceBuf_[i] = NULL;
			}
			
			int nNumThread = EnvVal::getNumThreads();
                        //VEO_DEBUG_PRINT2("nNumThread", nNumThread, "n_veo_thread_count", n_veo_thread_count);
			if (pEpEpFunc == NULL) {
                        	// 環境変数 OMP_NUM_THREADSから1プロセスあたりのスレッド数を得る
                        	for (int ith=0; ith<nNumThread; ith++) {
					int ith2 = ith % n_veo_thread_count;
#if defined(VEO_PROC_EACH_VECORE)
					pForceBuf_[ith] = pInnerVeoFuncEpEp_[ith2]->createVeBuffer(0);
#else
					pForceBuf_[ith] = pInnerVeoFuncEpEp_[ith2 / N_VE_NUM_CORES]->createVeBuffer(ith2 % N_VE_NUM_CORES);
#endif
				}

#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
				//p_omp_lock_ = new  omp_lock_t[n_veo_thread_count];
				p_omp_lock_ = (omp_lock_t *)malloc(sizeof(omp_lock_t) * n_veo_thread_count);
                        	for (int i=0; i<n_veo_thread_count; i++) {
                               		omp_init_lock(&p_omp_lock_[i]);
				}
#endif
			} else {
				bSuperParticleType_ = true;
                        	for (int i=0; i<N_VEO_MAX_THREAD; i++) {
					pForceBuf_[i] = pEpEpFunc->pForceBuf_[i];
				}

#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
				p_omp_lock_ = pEpEpFunc->p_omp_lock_;
#endif
			}

                        for (int ith=0; ith<nNumThread; ith++) {
				int ith2 = ith % n_veo_thread_count;
#if defined(VEO_PROC_EACH_VECORE)
				pEpjForForceBuf_[ith] = pInnerVeoFuncEpEp_[ith2]->createVeBuffer(0);
#else
				pEpjForForceBuf_[ith] = pInnerVeoFuncEpEp_[ith2 / N_VE_NUM_CORES]->createVeBuffer(ith2 % N_VE_NUM_CORES);
#endif
			}

#if defined(VEO_USE_PROFILE)
			prepare_profile_.symNum_ = 0;
			//strcpy(prepare_profile_.funcName_, psym_prepare_name);
			strcpy(prepare_profile_.funcName_, "prepare_");
			strcat(prepare_profile_.funcName_, psym_name);

			profile_.symNum_ = symNum;
			strcpy(profile_.funcName_, psym_name);
#endif
		}

		// デストラクタ
		~VeoFuncEpEp() {
#if defined(VEO_USE_PROFILE)
			prepare_profile_.output();
			profile_.output();
#endif
			if(!bSuperParticleType_) {
				for (int i=0; i<N_VEO_MAX_THREAD; i++) {
					if (pForceBuf_[i] != NULL) delete pForceBuf_[i];
				}
			}
			for (int i=0; i<N_VEO_MAX_THREAD; i++) {
				if (pEpjForForceBuf_[i] != NULL) delete pEpjForForceBuf_[i];
			}
			for (int i=0; i<N_VEO_MAX_VE; i++) {
				if (pInnerVeoFuncEpEp_[i] != NULL) {
					delete pInnerVeoFuncEpEp_[i];
				}
			}

			if (!bSuperParticleType_) {
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
                        	for (int i=0; i<n_veo_thread_count; i++) {
                               		omp_destroy_lock(&p_omp_lock_[i]);
				}

				//delete [] p_omp_lock_;
				free(p_omp_lock_);
#endif
			}
		}

		void setLongType() {
			bLongType_ = true;
		}

		void downloadSortedEpEp(const void* pSortedEpiBuf, size_t epi_size, const void* pSortedEpjBuf, size_t epj_size) {
#if defined(VEO_PROC_EACH_VECORE)
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
                        #pragma omp parallel
                        {
                                int ith = omp_get_thread_num();
                                //if (ith >= 0 && ith < n_veo_thread_count && pInnerVeoFuncEpEp_[ith] != NULL) {
                                if (ith < n_veo_thread_count) {
                                        pInnerVeoFuncEpEp_[ith]->downloadSortedEpEp(pSortedEpiBuf, epi_size, pSortedEpjBuf, epj_size);
                                }

                                #pragma omp barrier
                        }
#else
                        pInnerVeoFuncEpEp_[0]->downloadSortedEpEp(pSortedEpiBuf, epi_size, pSortedEpjBuf, epj_size);
#endif

#else // not VEO_PROC_EACH_VECORE
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
			#pragma omp parallel
			{
				int ith = omp_get_thread_num();
				if (ith < n_veo_thread_count && ith % N_VE_NUM_CORES == 0) {
					pInnerVeoFuncEpEp_[ith / N_VE_NUM_CORES]->downloadSortedEpEp(pSortedEpiBuf, epi_size, pSortedEpjBuf, epj_size);
				}
				
				#pragma omp barrier
			}
#else
			pInnerVeoFuncEpEp_[0]->downloadSortedEpEp(pSortedEpiBuf, epi_size, pSortedEpjBuf, epj_size);
#endif
#endif
		}

		void downloadSortedSp(const void* pSortedEpjBuf, size_t epj_size) {
#if defined(VEO_PROC_EACH_VECORE)
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
                        #pragma omp parallel
                        {
                                int ith = omp_get_thread_num();
                                //if (ith >= 0 && ith < n_veo_thread_count && pInnerVeoFuncEpEp_[ith] != NULL) {
                                if (ith < n_veo_thread_count) {
                                        pInnerVeoFuncEpEp_[ith]->downloadSortedSp(pSortedEpjBuf, epj_size);
                                }

                                #pragma omp barrier
                        }
#else
                        pInnerVeoFuncEpEp_[0]->downloadSortedSp(pSortedEpjBuf, epj_size);
#endif

#else // not VEO_PROC_EACH_VECORE
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL) && defined(_OPENMP)
			#pragma omp parallel
			{
				int ith = omp_get_thread_num();
				if (ith < n_veo_thread_count && ith % N_VE_NUM_CORES == 0) {
					pInnerVeoFuncEpEp_[ith / N_VE_NUM_CORES]->downloadSortedSp(pSortedEpjBuf, epj_size);
				}
				
				#pragma omp barrier
			}
#else
			pInnerVeoFuncEpEp_[0]->downloadSortedSp(pSortedEpjBuf, epj_size);
#endif
#endif
		}

// aaaaa

#if defined(PREPARE_EPJ_IN_VE)
		int callFuncEpEp(int ith, size_t offset_epi_size, int n_epi,
				const ParticleSimulator::S32* padr_epj_force_buf, size_t epj_size, int n_epj,
								void *pForceBuf, size_t force_size, const double *params = NULL, int params_count = 0) {
			int ret = -1;
			int ith0 = ith;
			ith = ith % n_veo_thread_count;
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
			omp_set_lock(&p_omp_lock_[ith]);
#endif
#if defined(VEO_PROC_EACH_VECORE)
			ret = pInnerVeoFuncEpEp_[ith]->callFuncEpEp(0, offset_epi_size, n_epi,
								padr_epj_force_buf, epj_size, n_epj,
								pEpjForForceBuf_[ith0],
								pForceBuf_[ith0], pForceBuf, force_size, params, params_count, bLongType_);
#else
			ret = pInnerVeoFuncEpEp_[ith / N_VE_NUM_CORES]->callFuncEpEp(ith % N_VE_NUM_CORES, offset_epi_size, n_epi,
								padr_epj_force_buf, epj_size, n_epj,
								pEpjForForceBuf_[ith0],
								pForceBuf_[ith0], pForceBuf, force_size, params, params_count, bLongType_);
#endif
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
			omp_unset_lock(&p_omp_lock_[ith]);
#endif

			return ret;
		}

#else	// not PREPARE_EPJ_IN_VE
                int callFuncEpEp(int ith, size_t offset_epi_size, int n_epi,
                                const void * epj_for_force_buf, size_t epj_size, int n_epj,
                                                                void *pForceBuf, size_t force_size, const double *params = NULL, int params_count = 0) {
                        int ret = -1;
                        int ith0 = ith;
                        ith = ith % n_veo_thread_count;
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
                        omp_set_lock(&p_omp_lock_[ith]);
#endif
#if defined(VEO_PROC_EACH_VECORE)
                        ret = pInnerVeoFuncEpEp_[ith]->callFuncEpEp(0, offset_epi_size, n_epi,
                                                                epj_for_force_buf, epj_size, n_epj,
                                                                pEpjForForceBuf_[ith0],
								pForceBuf_[ith0], pForceBuf, force_size, params, params_count, bLongType_);
#else
                        ret = pInnerVeoFuncEpEp_[ith / N_VE_NUM_CORES]->callFuncEpEp(ith % N_VE_NUM_CORES, offset_epi_size, n_epi,
                                                                epj_for_force_buf, epj_size, n_epj,
                                                                pEpjForForceBuf_[ith0],
								pForceBuf_[ith0], pForceBuf, force_size, params, params_count, bLongType_);
#endif
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
                        omp_unset_lock(&p_omp_lock_[ith]);
#endif

                        return ret;
                }
#endif	// end PREPARE_EPJ_IN_VE

#if defined(PREPARE_SPJ_IN_VE)
		int callFuncEpSp(int ith, size_t offset_epi_size, int n_epi,
				ParticleSimulator::S32* padr_spj_for_force_buf, size_t spj_size, int n_spj,
								void *pForceBuf, size_t force_size, const double *params = NULL, int params_count = 0) {
			int ret = -1;
			int ith0 = ith;
			ith = ith % n_veo_thread_count;
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
			omp_set_lock(&p_omp_lock_[ith]);
#endif
#if defined(VEO_PROC_EACH_VECORE)
			ret = pInnerVeoFuncEpEp_[ith]->callFuncEpSp(0, offset_epi_size, n_epi,
								padr_spj_for_force_buf, spj_size, n_spj,
                                                                pEpjForForceBuf_[ith0],
								pForceBuf_[ith0], pForceBuf, force_size, params, params_count);
#else
			ret = pInnerVeoFuncEpEp_[ith / N_VE_NUM_CORES]->callFuncEpSp(ith % N_VE_NUM_CORES, offset_epi_size, n_epi,
								padr_spj_for_force_buf, spj_size, n_spj,
                                                                pEpjForForceBuf_[ith0],
								pForceBuf_[ith0], pForceBuf, force_size, params, params_count);
#endif
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
			omp_unset_lock(&p_omp_lock_[ith]);
#endif

			return ret;
		}
#else	// not PREPARE_SPJ_IN_VE
		int callFuncEpSp(int ith, size_t offset_epi_size, int n_epi,
				const void* spj_for_force_buf, size_t spj_size, int n_spj,
								void *pForceBuf, size_t force_size, const double *params = NULL, int params_count = 0) {
			int ret = -1;
			int ith0 = ith;
			ith = ith % n_veo_thread_count;
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
			omp_set_lock(&p_omp_lock_[ith]);
#endif
#if defined(VEO_PROC_EACH_VECORE)
			ret = pInnerVeoFuncEpEp_[ith]->callFuncEpSp(0, offset_epi_size, n_epi,
								spj_for_force_buf, spj_size, n_spj,
                                                                pEpjForForceBuf_[ith0],
								pForceBuf_[ith0], pForceBuf, force_size, params, params_count);
#else
			ret = pInnerVeoFuncEpEp_[ith / N_VE_NUM_CORES]->callFuncEpSp(ith % N_VE_NUM_CORES, offset_epi_size, n_epi,
								spj_for_force_buf, spj_size, n_spj,
                                                                pEpjForForceBuf_[ith0],
								pForceBuf_[ith0], pForceBuf, force_size, params, params_count);
#endif
#if defined(PARTICLE_SIMULATOR_THREAD_PARALLEL)
			omp_unset_lock(&p_omp_lock_[ith]);
#endif

			return ret;
		}
#endif	// end PREPARE_SPJ_IN_VE
	};


	// 複数プロセスで使用するVEの番号を調停するクラス
	class VeNumberFixer {
        	// 内部クラス
		// VEのコアの仮割り当てを行う
        	class VeCoreAssumption {
                	private:
                        	int iVeNumber_;
                        	int nFreeCoreCount_ = N_VE_NUM_CORES;
                	public:
                        VeCoreAssumption(int iVeNumber) {
                                	iVeNumber_ = iVeNumber;
                        }

                        int getVeNumber() const { return iVeNumber_; }

                        bool allocateCore(int iCoreCount) {
                                bool ret = false;

                                if (nFreeCoreCount_ >= iCoreCount) {
                                        nFreeCoreCount_ -= iCoreCount;
                                        ret = true;
                                }

                                return ret;
                        }
        	};

        private:
                int iNumThreads_;
                int iVeNumber_[N_VEO_MAX_THREAD];

        public:
                VeNumberFixer() {
                	int nProc = 1;
                	int nRank = 0;
			for(int i=0; i < N_VEO_MAX_THREAD; i++) iVeNumber_[i] = 0;

#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
                        MPI_Comm_rank(MPI_COMM_WORLD, &nRank);   // ランクの取得
                        MPI_Comm_size(MPI_COMM_WORLD, &nProc);   // 全プロセス数の取得
#endif

                        // 環境変数 OMP_NUM_THREADSから1プロセスあたりのスレッド数を得る
			int nNumThread = EnvVal::getNumThreads();

                        // 環境変数 VE_NODE_NUMBERから使用するVEノード番号の範囲を取得する
                        std::vector<int> arrayVeNumber2 = EnvVal::getIntRange("VE_NODE_NUMBER", 0, 7);
                        int nNumThread2 = arrayVeNumber2.size() * N_VE_NUM_CORES / nProc;
                        if (nNumThread > nNumThread2) {
                        	iNumThreads_ = nNumThread = nNumThread2;
                        }
                        else {
                        	iNumThreads_ = nNumThread;
                        }
                        //VEO_DEBUG_PRINT2("iNumThreads_", iNumThreads_, "nNumThread", nNumThread);

                        if (nRank == 0) {
				// VH側の1スレッドがVEの1プロセスに対応するようにする。
				// VH側がOpenMPの8スレッド/プロセスならば VE側は8プロセス/VEとする 
				//VEO_DEBUG_PRINT("nProc", nProc);
				//VEO_DEBUG_PRINT("nNumThread", nNumThread);

                                // 環境変数 VE_NODE_NUMBERから使用するVEノード番号の範囲を取得する
                                std::vector<int> arrayVeNumber = EnvVal::getIntRange("VE_NODE_NUMBER", 0, 7);

                                std::vector<VeCoreAssumption> arrayVeCoreAssumption;
                                for(auto itr=arrayVeNumber.begin(); itr != arrayVeNumber.end(); ++itr) {
                                        arrayVeCoreAssumption.push_back(VeCoreAssumption(*itr));
                                }

                                for(int i=0; i < nProc ; i++) {
					int ve_number[N_VEO_MAX_THREAD];

					for(int j=0; j < nNumThread; j++) {
                                        	for(auto itr=arrayVeCoreAssumption.begin(); itr != arrayVeCoreAssumption.end(); ++itr) {
                                                	if (itr->allocateCore(1)) {
                                                        	ve_number[j] = itr->getVeNumber();
                                                        	break;
                                                	}
						}
                                        }

                                        if (i == 0) {
						for(int j=0; j < nNumThread; j++) {
                                                	iVeNumber_[j] = ve_number[j];
						}
                                        }
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
                                        else {
						// Rank=1以上のプロセスへ送信
                                                MPI_Send(&ve_number[0], nNumThread, MPI_INT, i, 0, MPI_COMM_WORLD);
                                        }
#endif
                                }
                        }
#if defined(PARTICLE_SIMULATOR_MPI_PARALLEL)
                        else {
				// Rank=0のプロセスから受信
                                MPI_Recv(&iVeNumber_[0], nNumThread, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        }
#endif
			
			//VEO_DEBUG_PRINT2("nRank", nRank, "iVeNumber_", iVeNumber_);
		}

                int getNumThreads() const { return iNumThreads_; }

                int getVeNumber(int ith) const {
			if (ith >= 0 && ith < N_VEO_MAX_THREAD) {
                        	return iVeNumber_[ith];
			} else {
				return -1;
			}
                }
	};

	class VeoInterface {
	private:
		static VeoInterface* pInstance_ ;

#if defined(VEO_PROC_EACH_VECORE)
		VeoProcess *pVeoProcess_[N_VEO_MAX_THREAD];
#else
		VeoProcess *pVeoProcess_[N_VEO_MAX_VE];
#endif

#if defined(VEO_PROC_EACH_VECORE)
                VeoInterface(const char* pLibPath) {
                        for(int i=0; i<N_VEO_MAX_THREAD;i++) {
                                pVeoProcess_[i] = NULL;
                        }

                        VeNumberFixer   vfix;

                        for(int i=0; i<vfix.getNumThreads();i++) {
                                pVeoProcess_[i] = new VeoProcess(pLibPath, vfix.getVeNumber(i), 1);
                        }
		}
#else
		// コンストラクタ
		VeoInterface(const char* pLibPath) {
			int ve_context_count[N_VEO_MAX_VE];
			for(int i=0; i<N_VEO_MAX_VE;i++) {
				ve_context_count[i] = 0;
				pVeoProcess_[i] = NULL;
			}

			VEO_DEBUG_PRINT2("veo_api_version()", veo_api_version(), "veo_version_string()", veo_version_string());

			VeNumberFixer	vfix;

			for(int i=0; i<vfix.getNumThreads();i++) {
				ve_context_count[vfix.getVeNumber(i)] += 1;
			}

			int j = 0;
			for(int i=0; i<N_VEO_MAX_VE;i++) {
				if (ve_context_count[i] > 0) {
					pVeoProcess_[j] = new VeoProcess(pLibPath, i, ve_context_count[i]);
					j += 1;
				}
			}
		}
#endif

	public:
#if defined(VEO_PROC_EACH_VECORE)
		// デストラクタ
                ~VeoInterface() {
                        for(int i=0; i<N_VEO_MAX_THREAD;i++) {
                                if (pVeoProcess_[i] != NULL) delete pVeoProcess_[i];
                        }
                }
#else
		// デストラクタ
		~VeoInterface() {
			for(int i=0; i<N_VEO_MAX_VE;i++) {
				if (pVeoProcess_[i] != NULL) delete pVeoProcess_[i];
			}
		}
#endif

		// VEOffloadingインタフェース初期化
		static void initialize(const char* argv0, const char* libname) {
			// フルパスを得る
			char  argv0_fullpath[PATH_MAX];
			char  *argv0_res = realpath(argv0, argv0_fullpath);

			// ダイナミックリンクライブラリが実行モジュールと同じディレクトリにあるとする
			std::string     veo_libpath(dirname(argv0_res));
			veo_libpath += "/";
			veo_libpath += libname;

			VEO_DEBUG_PRINT2("PATH_MAX", PATH_MAX, "veo_libpath", veo_libpath);

			//std::cout << "Initialize interface for VE Offloading version " << veo_version_string() << " ..." << std::endl;
			//fprintf(stdout, "libname=[%s]\n", libname);
			//fprintf(stdout, "veo_libpath=[%s]\n", veo_libpath.c_str());
			VeoInterface::pInstance_ = new VeoInterface(veo_libpath.c_str());	
		} 

		static VeoInterface* getInstance() { return VeoInterface::pInstance_; }

		static void disposeInterface() {
			VEO_DEBUG_PRINT("pInstance_", pInstance_);
			delete VeoInterface::pInstance_;
			VeoInterface::pInstance_ = NULL;
			//std::cout << "Teraminate interface for VE Offloading version " << veo_version_string() << std::endl;
		}

#if defined(VEO_USE_PROFILE)
		VeoFuncEpEp* createFuncEpEp(int symNum, const char* psym_prepare_name, const char* psym_name, VeoFuncEpEp* pBaseFuncEpEp) {
			//VEO_DEBUG_PRINT("pInstance_", pInstance_);
			return new VeoFuncEpEp(pVeoProcess_, symNum, psym_prepare_name, psym_name, pBaseFuncEpEp);
		}
#else
		VeoFuncEpEp* createFuncEpEp(const char* psym_prepare_name, const char* psym_name, VeoFuncEpEp* pBaseFuncEpEp) {
			//VEO_DEBUG_PRINT("pInstance_", pInstance_);
			return new VeoFuncEpEp(pVeoProcess_, 0, psym_prepare_name, psym_name, pBaseFuncEpEp);
		}
#endif
	};

	template<class Tepi, class Tepj>
	class FunctorEpEp {
	private:
		VeoFuncEpEp *pFuncEpEp_;
	protected:
#if defined(VEO_USE_PROFILE)
		FunctorEpEp(int symNum, const char* pFuncPrepareEpjSymbol, const char* pFuncEpEpSymbol) {
			pFuncEpEp_ = VeoInterface::getInstance()->createFuncEpEp(symNum, pFuncPrepareEpjSymbol, pFuncEpEpSymbol, NULL);
		}
		FunctorEpEp(int symNum, const char* pFuncEpEpSymbol) {
			pFuncEpEp_ = VeoInterface::getInstance()->createFuncEpEp(symNum, "dummy", pFuncEpEpSymbol, NULL);
		}
#else
		FunctorEpEp(const char* pFuncPrepareEpjSymbol, const char* pFuncEpEpSymbol) {
			pFuncEpEp_ = VeoInterface::getInstance()->createFuncEpEp(pFuncPrepareEpjSymbol, pFuncEpEpSymbol, NULL);
		}
		FunctorEpEp(const char* pFuncEpEpSymbol) {
			pFuncEpEp_ = VeoInterface::getInstance()->createFuncEpEp("dummy", pFuncEpEpSymbol, NULL);
		}
#endif
		virtual ~FunctorEpEp() {
			//delete pFuncEpEp_;
		}

		double *params_ = NULL;

		int num_of_params_ = 0;

#if defined(PREPARE_EPJ_IN_VE)
		void perform(int ith, size_t offset_epi_size, int n_epi,
				const ParticleSimulator::S32* padr_epj_force_buf, size_t epj_size, int n_epj,
								void *pForceBuf, size_t force_size) {
			pFuncEpEp_->callFuncEpEp(ith, offset_epi_size, n_epi,
							padr_epj_force_buf, epj_size, n_epj,
							pForceBuf, force_size, params_, num_of_params_);
		}
#else
                void perform(int ith, size_t offset_epi_size, int n_epi,
                                const void* epj_for_force_buf, size_t epj_size, int n_epj,
                                                                void *pForceBuf, size_t force_size) {
                        pFuncEpEp_->callFuncEpEp(ith, offset_epi_size, n_epi,
                                                        epj_for_force_buf, epj_size, n_epj,
                                                        pForceBuf, force_size, params_, num_of_params_);
                }
#endif

		void transfer(const Tepi* const ep_i, const ParticleSimulator::S32 Nip,
                     const Tepj* const ep_j, const ParticleSimulator::S32 Njp) {
			pFuncEpEp_->downloadSortedEpEp((const void *)ep_i, sizeof(Tepi) * Nip, (const void *)ep_j, sizeof(Tepj) * Njp);
		} 

		void terminateFunctorEpEp() {
			delete pFuncEpEp_;
			if (params_ != NULL) {
				delete [] params_;
			}
		}
	public:
#if defined(PREPARE_EPJ_IN_VE)
		void operator() (int ith, size_t offset_epi_size, int n_epi,
			ParticleSimulator::S32 * padr_epj_for_force, size_t epj_size, int n_epj,
			void *pForceBuf, size_t force_size) {
			perform(ith, offset_epi_size, n_epi, padr_epj_for_force, epj_size, n_epj, pForceBuf, force_size);
		}
#else
		void operator() (int ith, size_t offset_epi_size, int n_epi,
			const void * epj_for_force_buf, size_t epj_size, int n_epj,
			void *pForceBuf, size_t force_size) {
			perform(ith, offset_epi_size, n_epi, epj_for_force_buf, epj_size, n_epj, pForceBuf, force_size);
		}
#endif

		void transferEpiEpj(const Tepi* const ep_i, const ParticleSimulator::S32 Nip,
			const Tepj* const ep_j, const ParticleSimulator::S32 Njp) {
			transfer(ep_i, Nip, ep_j, Njp);
		}

		void cleanup() {
			terminateFunctorEpEp();
		}

                VeoFuncEpEp* getFuncEpEpPtr() const {
                        return pFuncEpEp_;
                }

		void setLongType() {
			pFuncEpEp_->setLongType();
		}
	};


	template<class Tepi, class Tepj, class Tspj>
	class FunctorEpSp {
	private:
		VeoFuncEpEp *pFuncEpEp_;
	protected:
#if defined(PREPARE_SPJ_IN_VE)
#if defined(VEO_USE_PROFILE)
		FunctorEpSp(int symNum, const char* pFuncPrepareSpjSymbol, const char* pFuncEpEpSymbol, FunctorEpEp<Tepi,Tepj>& functorEpEp) {
			pFuncEpEp_ = VeoInterface::getInstance()->createFuncEpEp(symNum, pFuncPrepareSpjSymbol, pFuncEpEpSymbol, functorEpEp.getFuncEpEpPtr());
			functorEpEp.setLongType();
		}
#else
		FunctorEpSp(const char* pFuncPrepareSpjSymbol, const char* pFuncEpEpSymbol, FunctorEpEp<Tepi,Tepj>& functorEpEp ) {
			pFuncEpEp_ = VeoInterface::getInstance()->createFuncEpEp(pFuncPrepareSpjSymbol, pFuncEpEpSymbol, functorEpEp.getFuncEpEpPtr());
			functorEpEp.setLongType();
		}
#endif
#else
#if defined(VEO_USE_PROFILE)
		FunctorEpSp(int symNum, const char* pFuncEpEpSymbol, FunctorEpEp<Tepi,Tepj>& functorEpEp) {
			pFuncEpEp_ = VeoInterface::getInstance()->createFuncEpEp(symNum, "dummy", pFuncEpEpSymbol, functorEpEp.getFuncEpEpPtr());
			functorEpEp.setLongType();
		}
		FunctorEpSp(int symNum, const char* pFuncPrepareSpjSymbol, const char* pFuncEpEpSymbol, FunctorEpEp<Tepi,Tepj>& functorEpEp) {
			pFuncEpEp_ = VeoInterface::getInstance()->createFuncEpEp(symNum, "dummy", pFuncEpEpSymbol, functorEpEp.getFuncEpEpPtr());
			functorEpEp.setLongType();
		}
#else
		FunctorEpSp(const char* pFuncEpEpSymbol, FunctorEpEp<Tepi,Tepj>& functorEpEp ) {
			pFuncEpEp_ = VeoInterface::getInstance()->createFuncEpEp("dummy", pFuncEpEpSymbol, functorEpEp.getFuncEpEpPtr());
			functorEpEp.setLongType();
		}
		FunctorEpSp(const char* pFuncPrepareSpjSymbol, const char* pFuncEpEpSymbol, FunctorEpEp<Tepi,Tepj>& functorEpEp ) {
			pFuncEpEp_ = VeoInterface::getInstance()->createFuncEpEp("dummy", pFuncEpEpSymbol, functorEpEp.getFuncEpEpPtr());
			functorEpEp.setLongType();
		}
#endif
#endif
		virtual ~FunctorEpSp() {
			//delete pFuncEpEp_;
		}

		double *params_ = NULL;

		int num_of_params_ = 0;

#if defined(PREPARE_SPJ_IN_VE)
		void performSp(int ith, size_t offset_epi_size, int n_epi,
				ParticleSimulator::S32* padr_spj_for_force_buf, size_t spj_size, int n_spj,
								void *pForceBuf, size_t force_size, const double *params = NULL, int params_count = 0) {
			pFuncEpEp_->callFuncEpSp(ith, offset_epi_size, n_epi,
							padr_spj_for_force_buf, spj_size, n_spj,
							pForceBuf, force_size, params_, num_of_params_);
		}
#else
		void performSp(int ith, size_t offset_epi_size, int n_epi,
				const void* spj_for_force_buf, size_t spj_size, int n_spj,
								void *pForceBuf, size_t force_size, const double *params = NULL, int params_count = 0) {
			pFuncEpEp_->callFuncEpSp(ith, offset_epi_size, n_epi,
							spj_for_force_buf, spj_size, n_spj,
							pForceBuf, force_size, params_, num_of_params_);
		}
#endif

		void transferSp(const Tspj* const sp_j, const ParticleSimulator::S32 Njp) {
			pFuncEpEp_->downloadSortedSp((const void *)sp_j, sizeof(Tspj) * Njp);
		} 

		void terminateFunctorEpEp() {
			delete pFuncEpEp_;
		}

	public:
#if defined(PREPARE_SPJ_IN_VE)
		void operator() (int ith, size_t offset_epi_size, int n_epi,
			ParticleSimulator::S32* padr_spj_for_force_buf, size_t spj_size, int n_spj,
			void *pForceBuf, size_t force_size) {
			performSp(ith, offset_epi_size, n_epi, padr_spj_for_force_buf, spj_size, n_spj, pForceBuf, force_size);
		}
#else
		void operator() (int ith, size_t offset_epi_size, int n_epi,
			const void* spj_for_force_buf, size_t spj_size, int n_spj,
			void *pForceBuf, size_t force_size) {
			performSp(ith, offset_epi_size, n_epi, spj_for_force_buf, spj_size, n_spj, pForceBuf, force_size);
		}
#endif

		void transferSpj(const Tspj* const ep_j, const ParticleSimulator::S32 Njp) {
			transferSp(ep_j, Njp);
		}

		void cleanup() {
			terminateFunctorEpEp();
		}
	};

#if !defined(_SX_FORTRAN_VEO)
	VeoInterface* VeoInterface::pInstance_ = NULL;
#endif

}; // end namespace VEOffloading

namespace VEO = VEOffloading;

//#define DECLARE_VEOFFLOAD_INTERFACE VEOffloading::VeoInterface *VEOffloading::VeoInterface::pInstance_=NULL
