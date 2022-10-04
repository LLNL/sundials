


struct ARKodeMethodMem {
  /* Further method specification */
  void* impl;

  int q;                      /* method order                                       */
  int p;                      /* embedding order (0 = no embedding)                 */
  int stages;                 /* number of stages                                   */
};

typedef struct ARKodeMethodMem* ARKodeMethod;
