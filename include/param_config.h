#ifndef PARAM_CONFIG_H
#define PARAM_CONFIG_H

#include <stdbool.h> // For bool type
#include <stddef.h>  // For size_t

// Enum for parameter types
typedef enum {
    PARAM_TYPE_INT,
    PARAM_TYPE_FLOAT,
    PARAM_TYPE_STRING, // For char arrays
    PARAM_TYPE_BOOL    // For bool variables
} ParamType;

// Structure to hold parameter metadata
typedef struct {
    const char *name;
    ParamType type;
    void *ptr; // Pointer to the actual variable
    size_t size; // Size for string types (MAXLENLINE, etc.)
    bool is_old_version; // Flag to indicate if this is an old/deprecated parameter name
} ParameterEntry;

#endif // PARAM_CONFIG_H