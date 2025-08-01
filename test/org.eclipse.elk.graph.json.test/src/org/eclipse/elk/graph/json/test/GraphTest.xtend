/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.graph.json.test

import com.google.gson.JsonIOException
import org.eclipse.elk.alg.test.PlainJavaInitialization
import org.eclipse.elk.core.options.CoreOptions
import org.eclipse.elk.core.options.Direction
import org.eclipse.elk.graph.json.ElkGraphJson
import org.eclipse.elk.graph.json.JsonImportException
import org.junit.BeforeClass
import org.junit.Test

import static org.junit.Assert.*

/**
 */
class GraphTest {
    
    @BeforeClass
    static def void init() {
        PlainJavaInitialization.initializePlainJavaLayout
    }
    
    @Test 
    def void graphMustBeObjectTest() {
        ElkGraphJson.forGraph("{\"id\":1}").toElk
    }
    
    @Test(expected = JsonImportException) 
    def void graphMustBeObjectFailArrayTest() {
        ElkGraphJson.forGraph("[]").toElk
    }
    
    @Test
    def void smallGraphTest() {
        val graph = '''
        {
          "id": "root",
          "properties": {
            "elk.direction": "DOWN"
          },
          "children": [{"id": "n1", "width": 40, "height": 40},
                      {"id": "n2", "width": 40, "height": 40}],
          "edges": [{"id": "e1", "source": "n1", "target": "n2"}]
        }
        '''
        
        val root = ElkGraphJson.forGraph(graph).toElk
        
        assertTrue(root.children.size === 2)
        assertTrue(root.containedEdges.size === 1)
        
        assertTrue(root.getProperty(CoreOptions.DIRECTION) == Direction.DOWN)
    }
    
    val sloppyJsonGraph = '''
        {
          // the root node
          id: "root",
          /* Now the graph */ 
          "children": [ {"id": "c"}; {"id": "c1"} ],
          'edges': [] // Endline comment
        }
    '''
    
    @Test(expected = JsonIOException)
    def void rejectSloppyJsonIfNotLenient() {
         ElkGraphJson.forGraph(sloppyJsonGraph)
                     .lenient(false)
                     .toElk
    }
    
    @Test
    def void acceptSloppyJsonIfLenient() {
        val root = ElkGraphJson.forGraph(sloppyJsonGraph).toElk
        
        assertEquals("root", root.identifier)
        assertEquals(2, root.children.size)
        assertEquals("c", root.children.head.identifier)
        assertEquals("c1", root.children.lastOrNull.identifier)
    }
}