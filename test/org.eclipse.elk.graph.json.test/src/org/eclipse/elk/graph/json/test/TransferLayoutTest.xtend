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

import com.google.gson.JsonObject
import com.google.gson.JsonParser
import org.eclipse.elk.core.RecursiveGraphLayoutEngine
import org.eclipse.elk.core.util.BasicProgressMonitor
import org.eclipse.elk.core.util.Maybe
import org.eclipse.elk.graph.json.ElkGraphJson
import org.eclipse.elk.graph.json.JsonImportException
import org.eclipse.elk.graph.json.JsonImporter
import org.eclipse.elk.graph.util.ElkGraphUtil
import org.junit.Test

import static org.junit.Assert.*

/**
 */
class TransferLayoutTest {
    
    val graph = '''
    {
      "id": "root",
      "properties": {
        "algorithm": "random",
        "direction": "DOWN"
      },
      "children": [{"id": "n1", "width": 40, "height": 40},
                  {"id": "n2", "width": 40, "height": 40}],
      "edges": [{"id": "e1", "source": "n1", "target": "n2"}]
    }
    '''
        
    @Test
    def void transferLayoutTest() {
        val mby = new Maybe<JsonImporter>
        val root = ElkGraphJson.forGraph(graph)
                               .rememberImporter(mby)
                               .toElk
        
        mby.get.transferLayout(root)
    }
    
    @Test(expected = JsonImportException)
    def void transferLayoutFailTest() {
        val mby = new Maybe<JsonImporter>
        val root = ElkGraphJson.forGraph(graph).rememberImporter(mby).toElk
        
        // modify graph
        root.children += ElkGraphUtil.createNode(root)
         
        // this should throw an exception now
        mby.get.transferLayout(root)
    }
    
    @Test
    def void transferLayoutEdgeSectionExistsTest() {
        val jsonGraph = JsonParser.parseString(graph).asJsonObject
        
        val mby = new Maybe<JsonImporter>
        val root = ElkGraphJson.forGraph(jsonGraph)
                               .rememberImporter(mby)
                               .toElk
        
        new RecursiveGraphLayoutEngine().layout(root, new BasicProgressMonitor)
        
        mby.get.transferLayout(root)

        jsonGraph.get("children").asJsonArray.filter(JsonObject).forEach[ c |
            c.checkCoordinate("x")
            c.checkCoordinate("y")
        ]
        
        jsonGraph.get("edges").asJsonArray.filter(JsonObject).forEach[ e |
            val sections = e.get("sections")
            assertTrue("Every edge must have at least one section.", sections !== null)
            sections.asJsonArray.filter(JsonObject).forEach[checkSection]
            e.checkContainer();
        ]
        
    }
    
    private def checkCoordinate(JsonObject o, String c) {
        assertTrue(c + " coordinate must exist", o.has(c))
        val coor = o.get(c)
        assertTrue(coor.isJsonPrimitive)
        assertTrue(coor.asJsonPrimitive.asDouble > 0)
    }
    
    private def checkSection(JsonObject o) {
        assertTrue(o.has("id"))
        assertTrue(o.has("startPoint"))
        assertTrue(o.has("endPoint"))
    }

    private def checkContainer(JsonObject o) {
        assertTrue(o.has("container"))
    }
}