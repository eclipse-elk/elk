/**
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 */
package org.eclipse.elk.graph.json

import com.google.gson.JsonArray
import com.google.gson.JsonElement
import com.google.gson.JsonObject
import com.google.gson.JsonPrimitive

/**
 * The library dependent part of the {@link JsonImporter} 
 * using Google's gson library.
 * 
 * <b>IMPORTANT</b>
 * Whenever methods are added or removed from this class, 
 * elkjs's corresponding JsonAdapter must be update accordingly!
 */
final class JsonAdapter {

    // - - - - - - - - - - - - - - - - - - -
    // Getting the id of an element
    // - - - - - - - - - - - - - - - - - - -
    
    def Object getId(JsonObject o) {
        if (!o.has("id")) {
            throw formatError("Every element must have an id.")            
        }
        o.get("id").asId
    }
    
    def Object asId(JsonElement id) {
        switch id {
            JsonPrimitive case id.isString: return id.asString
            JsonPrimitive case id.isNumber && isInt(id.asDouble): return id.asDouble.intValue
            default: throw formatError("Id must be a string or an integer: '" + id + "'.")
        }   
    }
    
    def String getIdSave(JsonObject o) {
        o.optString("id")
    }
    
    // - - - - - - - - - - - - - - - - - - -
    // Errors
    // - - - - - - - - - - - - - - - - - - -
    
    def implementationError() {
        new Error("Severe implementation error in the Json to ElkGraph importer.")
    }
    
    def formatError(String msg) {
        new JsonImportException(msg)
    }

    // - - - - - - - - - - - - - - - - - - -    
    // Type emulation
    // - - - - - - - - - - - - - - - - - - -
    
    def JsonObject toJsonObject(Object o) {
        return o as JsonObject
    }
    
    def JsonArray toJsonArray(Object o) {
        return o as JsonArray
    }
    
    // - - - - - - - - - - - - - - - - - - -
    // Getting values of the json elements
    // - - - - - - - - - - - - - - - - - - -
    
    /**
     * If {@value e} is a primitve type (string, number, boolean) 
     * it is guaranteed to be converted to a String
     * (in the gson implementation of 'asString').
     */
    def String stringVal(JsonElement e) {
        e.asString
    }
    
    def String optString(JsonObject o, String element) {
        o.get(element)?.stringVal
    }
    
    def String optString(JsonArray arr, int i) {
        arr.get(i)?.stringVal
    }
    
    def String optString(Object o, String element) {
        o?.toJsonObject.optString(element)
    }

    def Double optDouble(JsonObject o, String element) {
        if (o.has(element))
            o.get(element).asDouble
        else 
            null
    }

    def JsonArray optJSONArray(JsonObject o, String element) {
        o.get(element)?.asJsonArray
    }

    def JsonObject optJSONObject(JsonObject o, String element) {
        o.get(element)?.asJsonObject
    }

    def JsonObject optJSONObject(JsonArray arr, int i) {
        arr.get(i)?.asJsonObject
    }
    
    // - - - - - - - - - - - - - - - - - - -
    
    def int sizeJsonArr(JsonArray o) {
        o.size
    }
    
    def boolean hasJsonObj(JsonObject o, String element) {
        o.has(element)
    }
    
    def Iterable<String> keysJsonObj(JsonObject o) {
        o.entrySet.map[e | e.key ]
    }
    
    def JsonElement getJsonObj(JsonObject o, String element) {
        o.get(element)
    }
    
    def JsonElement getJsonArr(JsonArray arr, int i) {
        arr.get(i)
    }
    
    // - - - - - - - - - - - - - - - - - - -
    // Create new json elements
    // - - - - - - - - - - - - - - - - - - -
    
    def JsonObject newJsonObject() {
        return new JsonObject
    }
    
    def JsonArray newJsonArray() {
        return new JsonArray
    }
    
    def JsonPrimitive toJson(String s) {
        return new JsonPrimitive(s)
    }
    
    def JsonPrimitive toJson(Number n) {
        return new JsonPrimitive(n)
    }
    
    // - - - - - - - - - - - - - - - - - - -
    // Add members to a json object
    // - - - - - - - - - - - - - - - - - - -
    
    def void addJsonObj(JsonObject o, String element, Number n) {
        o.addProperty(element, n)
    }
    
    def void addJsonObj(JsonObject o, String element, boolean b) {
        o.addProperty(element, b)
    }
    
    def void addJsonObj(JsonObject o, String element, String s) {
        o.addProperty(element, s)
    }
    
    def void addJsonObj(JsonObject o, String element, JsonElement je) {
        o.add(element, je)
    }
    
    def void addJsonObj(JsonObject o, String element, Object obj) {
        switch obj {
            String: o.addJsonObj(element, obj)
            Boolean: o.addJsonObj(element, obj)
            Number: o.addJsonObj(element, obj)
            default: throw implementationError
        }
    }
    
    // - - - - - - - - - - - - - - - - - - -
    // Remove members from a json object
    // - - - - - - - - - - - - - - - - - - -
    
    def void removeJsonObj(JsonObject o, String element) {
        o.remove(element)
    }
    
    // - - - - - - - - - - - - - - - - - - -
    // Add to a json array
    // - - - - - - - - - - - - - - - - - - -
        
    def void addJsonArr(JsonArray arr, JsonElement je) {
        arr.add(je)
    }

    def void addJsonArr(JsonArray arr, Object o) {
        switch o {
            String: arr.addJsonArr(o.toJson)
            Number: arr.addJsonArr(o.toJson)
            default: throw implementationError
        }
    }
    
    // - - - - - - - - - - - - - - - - - - -
    // Convenience
    // - - - - - - - - - - - - - - - - - - -
    
    def isInt(double d) {
       return d % 1 === 0
    }
}
