/*******************************************************************************
 * Copyright (c) 2024 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.graph.json.test

import com.google.gson.JsonParser
import org.eclipse.elk.alg.test.PlainJavaInitialization
import org.eclipse.elk.core.RecursiveGraphLayoutEngine
import org.eclipse.elk.core.util.BasicProgressMonitor
import org.eclipse.elk.core.util.Maybe
import org.eclipse.elk.graph.json.ElkGraphJson
import org.eclipse.elk.graph.json.JsonImporter
import org.junit.BeforeClass
import org.junit.Test

import static org.junit.Assert.*

/**
 */
class EdgeCoordsTest {
        
    @BeforeClass
    static def void init() {
        PlainJavaInitialization.initializePlainJavaLayout
    }
    val graph = '''
    {
        "id": "root",
        "properties": {
            "algorithm": "layered",
            "org.eclipse.elk.hierarchyHandling": "INCLUDE_CHILDREN"
        },
        "children": [
            { "id": "A",
                "children": [
                    { "id": "x", "width": 50, "height": 90 },
                    { "id": "B",
                        "labels": [ { "text": "B", "width": 10, "height": 12 } ],
                        "ports": [
                            { "id": "p", "width": 10, "height": 10,
                                "labels": [ { "text": "p", "width": 10, "height": 12 } ]
                            }
                        ],
                        "children": [
                            { "id": "y", "width": 50, "height": 90 },
                            { "id": "z", "width": 50, "height": 90 }
                        ],
                        "edges": [
                            { "id": "e1", "sources": [ "y" ], "targets": [ "z" ] },
                            { "id": "e2", "sources": [ "x" ], "targets": [ "z" ],
                                "labels": [ { "text": "e2", "width": 20, "height": 12 } ]
                            },
                            { "id": "e3", "sources": [ "x" ], "targets": [ "p" ] },
                            { "id": "e4", "sources": [ "p" ], "targets": [ "y" ] }
                        ]
                    }
                ]
            }
        ]
    }
    '''
            
    @Test
    def void edgeCoordsTest() {
        
        val cases = #[
            #['PARENT', 'CONTAINER', outputPC],
            #['PARENT', 'PARENT', outputPP],
            #['PARENT', 'ROOT', outputPR],
            #['ROOT', 'CONTAINER', outputRC],
            #['ROOT', 'PARENT', outputRP],
            #['ROOT', 'ROOT', outputRR]
        ]
        
        for (c : cases) {
            val scm = c.get(0)
            val ecm = c.get(1)
            val expectedOutputString = c.get(2)
            
            val jsonGraph = JsonParser.parseString(graph).asJsonObject
            
            val obj = jsonGraph.get("properties").asJsonObject
            obj.addProperty(
                "org.eclipse.elk.json.shapeCoords", scm
            )
            obj.addProperty(
                "org.eclipse.elk.json.edgeCoords", ecm
            )
            
            val mby = new Maybe<JsonImporter>
            val root = ElkGraphJson.forGraph(jsonGraph)
                                   .rememberImporter(mby)
                                   .toElk
            
            new RecursiveGraphLayoutEngine().layout(root, new BasicProgressMonitor)
            
            mby.get.transferLayout(root)
                           
            val computedOutput = jsonGraph.asJsonObject
            
            val expectedOutput = JsonParser.parseString(expectedOutputString).asJsonObject
            
            assertEquals(expectedOutput, computedOutput)
        }

    }
    
    val outputPC = '''
    {
     "id": "root",
     "properties": {
      "algorithm": "layered",
      "org.eclipse.elk.hierarchyHandling": "INCLUDE_CHILDREN",
      "org.eclipse.elk.json.shapeCoords": "PARENT",
      "org.eclipse.elk.json.edgeCoords": "CONTAINER"
     },
     "children": [
      {
       "id": "A",
       "children": [
        {
         "id": "x",
         "width": 50,
         "height": 90,
         "x": 12,
         "y": 39
        },
        {
         "id": "B",
         "labels": [
          {
           "text": "B",
           "width": 10,
           "height": 12,
           "x": 0,
           "y": 0
          }
         ],
         "ports": [
          {
           "id": "p",
           "width": 10,
           "height": 10,
           "labels": [
            {
             "text": "p",
             "width": 10,
             "height": 12,
             "x": -11,
             "y": -13
            }
           ],
           "x": -10,
           "y": 52
          }
         ],
         "children": [
          {
           "id": "y",
           "width": 50,
           "height": 90,
           "x": 12,
           "y": 12
          },
          {
           "id": "z",
           "width": 50,
           "height": 90,
           "x": 82,
           "y": 27
          }
         ],
         "edges": [
          {
           "id": "e1",
           "sources": [
            "y"
           ],
           "targets": [
            "z"
           ],
           "sections": [
            {
             "id": "e1_s0",
             "startPoint": {
              "x": 62,
              "y": 57
             },
             "endPoint": {
              "x": 82,
              "y": 57
             },
             "incomingShape": "y",
             "outgoingShape": "z"
            }
           ],
           "container": "B"
          },
          {
           "id": "e2",
           "sources": [
            "x"
           ],
           "targets": [
            "z"
           ],
           "labels": [
            {
             "text": "e2",
             "width": 20,
             "height": 12,
             "x": 82,
             "y": 102
            }
           ],
           "sections": [
            {
             "id": "e2_s0",
             "startPoint": {
              "x": 62,
              "y": 99
             },
             "endPoint": {
              "x": 225,
              "y": 99
             },
             "bendPoints": [
              {
               "x": 112,
               "y": 99
              },
              {
               "x": 112,
               "y": 124
              },
              {
               "x": 215,
               "y": 124
              },
              {
               "x": 215,
               "y": 99
              }
             ],
             "incomingShape": "x",
             "outgoingShape": "z"
            }
           ],
           "container": "A"
          },
          {
           "id": "e3",
           "sources": [
            "x"
           ],
           "targets": [
            "p"
           ],
           "sections": [
            {
             "id": "e3_s0",
             "startPoint": {
              "x": 62,
              "y": 69
             },
             "endPoint": {
              "x": 133,
              "y": 69
             },
             "incomingShape": "x",
             "outgoingShape": "p"
            }
           ],
           "container": "A"
          },
          {
           "id": "e4",
           "sources": [
            "p"
           ],
           "targets": [
            "y"
           ],
           "sections": [
            {
             "id": "e4_s0",
             "startPoint": {
              "x": -10,
              "y": 57
             },
             "endPoint": {
              "x": 12,
              "y": 57
             },
             "incomingShape": "p",
             "outgoingShape": "y"
            }
           ],
           "container": "B"
          }
         ],
         "x": 143,
         "y": 12,
         "width": 144,
         "height": 129
        }
       ],
       "x": 12,
       "y": 12,
       "width": 299,
       "height": 153
      }
     ],
     "x": 0,
     "y": 0,
     "width": 323,
     "height": 177
    }
    '''
    
    val outputPP = '''
    {
     "id": "root",
     "properties": {
      "algorithm": "layered",
      "org.eclipse.elk.hierarchyHandling": "INCLUDE_CHILDREN",
      "org.eclipse.elk.json.shapeCoords": "PARENT",
      "org.eclipse.elk.json.edgeCoords": "PARENT"
     },
     "children": [
      {
       "id": "A",
       "children": [
        {
         "id": "x",
         "width": 50,
         "height": 90,
         "x": 12,
         "y": 39
        },
        {
         "id": "B",
         "labels": [
          {
           "text": "B",
           "width": 10,
           "height": 12,
           "x": 0,
           "y": 0
          }
         ],
         "ports": [
          {
           "id": "p",
           "width": 10,
           "height": 10,
           "labels": [
            {
             "text": "p",
             "width": 10,
             "height": 12,
             "x": -11,
             "y": -13
            }
           ],
           "x": -10,
           "y": 52
          }
         ],
         "children": [
          {
           "id": "y",
           "width": 50,
           "height": 90,
           "x": 12,
           "y": 12
          },
          {
           "id": "z",
           "width": 50,
           "height": 90,
           "x": 82,
           "y": 27
          }
         ],
         "edges": [
          {
           "id": "e1",
           "sources": [
            "y"
           ],
           "targets": [
            "z"
           ],
           "sections": [
            {
             "id": "e1_s0",
             "startPoint": {
              "x": 62,
              "y": 57
             },
             "endPoint": {
              "x": 82,
              "y": 57
             },
             "incomingShape": "y",
             "outgoingShape": "z"
            }
           ]
          },
          {
           "id": "e2",
           "sources": [
            "x"
           ],
           "targets": [
            "z"
           ],
           "labels": [
            {
             "text": "e2",
             "width": 20,
             "height": 12,
             "x": -61,
             "y": 90
            }
           ],
           "sections": [
            {
             "id": "e2_s0",
             "startPoint": {
              "x": -81,
              "y": 87
             },
             "endPoint": {
              "x": 82,
              "y": 87
             },
             "bendPoints": [
              {
               "x": -31,
               "y": 87
              },
              {
               "x": -31,
               "y": 112
              },
              {
               "x": 72,
               "y": 112
              },
              {
               "x": 72,
               "y": 87
              }
             ],
             "incomingShape": "x",
             "outgoingShape": "z"
            }
           ]
          },
          {
           "id": "e3",
           "sources": [
            "x"
           ],
           "targets": [
            "p"
           ],
           "sections": [
            {
             "id": "e3_s0",
             "startPoint": {
              "x": -81,
              "y": 57
             },
             "endPoint": {
              "x": -10,
              "y": 57
             },
             "incomingShape": "x",
             "outgoingShape": "p"
            }
           ]
          },
          {
           "id": "e4",
           "sources": [
            "p"
           ],
           "targets": [
            "y"
           ],
           "sections": [
            {
             "id": "e4_s0",
             "startPoint": {
              "x": -10,
              "y": 57
             },
             "endPoint": {
              "x": 12,
              "y": 57
             },
             "incomingShape": "p",
             "outgoingShape": "y"
            }
           ]
          }
         ],
         "x": 143,
         "y": 12,
         "width": 144,
         "height": 129
        }
       ],
       "x": 12,
       "y": 12,
       "width": 299,
       "height": 153
      }
     ],
     "x": 0,
     "y": 0,
     "width": 323,
     "height": 177
    }
    '''
    
    val outputPR = '''
    {
     "id": "root",
     "properties": {
      "algorithm": "layered",
      "org.eclipse.elk.hierarchyHandling": "INCLUDE_CHILDREN",
      "org.eclipse.elk.json.shapeCoords": "PARENT",
      "org.eclipse.elk.json.edgeCoords": "ROOT"
     },
     "children": [
      {
       "id": "A",
       "children": [
        {
         "id": "x",
         "width": 50,
         "height": 90,
         "x": 12,
         "y": 39
        },
        {
         "id": "B",
         "labels": [
          {
           "text": "B",
           "width": 10,
           "height": 12,
           "x": 0,
           "y": 0
          }
         ],
         "ports": [
          {
           "id": "p",
           "width": 10,
           "height": 10,
           "labels": [
            {
             "text": "p",
             "width": 10,
             "height": 12,
             "x": -11,
             "y": -13
            }
           ],
           "x": -10,
           "y": 52
          }
         ],
         "children": [
          {
           "id": "y",
           "width": 50,
           "height": 90,
           "x": 12,
           "y": 12
          },
          {
           "id": "z",
           "width": 50,
           "height": 90,
           "x": 82,
           "y": 27
          }
         ],
         "edges": [
          {
           "id": "e1",
           "sources": [
            "y"
           ],
           "targets": [
            "z"
           ],
           "sections": [
            {
             "id": "e1_s0",
             "startPoint": {
              "x": 217,
              "y": 81
             },
             "endPoint": {
              "x": 237,
              "y": 81
             },
             "incomingShape": "y",
             "outgoingShape": "z"
            }
           ]
          },
          {
           "id": "e2",
           "sources": [
            "x"
           ],
           "targets": [
            "z"
           ],
           "labels": [
            {
             "text": "e2",
             "width": 20,
             "height": 12,
             "x": 94,
             "y": 114
            }
           ],
           "sections": [
            {
             "id": "e2_s0",
             "startPoint": {
              "x": 74,
              "y": 111
             },
             "endPoint": {
              "x": 237,
              "y": 111
             },
             "bendPoints": [
              {
               "x": 124,
               "y": 111
              },
              {
               "x": 124,
               "y": 136
              },
              {
               "x": 227,
               "y": 136
              },
              {
               "x": 227,
               "y": 111
              }
             ],
             "incomingShape": "x",
             "outgoingShape": "z"
            }
           ]
          },
          {
           "id": "e3",
           "sources": [
            "x"
           ],
           "targets": [
            "p"
           ],
           "sections": [
            {
             "id": "e3_s0",
             "startPoint": {
              "x": 74,
              "y": 81
             },
             "endPoint": {
              "x": 145,
              "y": 81
             },
             "incomingShape": "x",
             "outgoingShape": "p"
            }
           ]
          },
          {
           "id": "e4",
           "sources": [
            "p"
           ],
           "targets": [
            "y"
           ],
           "sections": [
            {
             "id": "e4_s0",
             "startPoint": {
              "x": 145,
              "y": 81
             },
             "endPoint": {
              "x": 167,
              "y": 81
             },
             "incomingShape": "p",
             "outgoingShape": "y"
            }
           ]
          }
         ],
         "x": 143,
         "y": 12,
         "width": 144,
         "height": 129
        }
       ],
       "x": 12,
       "y": 12,
       "width": 299,
       "height": 153
      }
     ],
     "x": 0,
     "y": 0,
     "width": 323,
     "height": 177
    }
    '''
    
    val outputRC = '''
    {
     "id": "root",
     "properties": {
      "algorithm": "layered",
      "org.eclipse.elk.hierarchyHandling": "INCLUDE_CHILDREN",
      "org.eclipse.elk.json.shapeCoords": "ROOT",
      "org.eclipse.elk.json.edgeCoords": "CONTAINER"
     },
     "children": [
      {
       "id": "A",
       "children": [
        {
         "id": "x",
         "width": 50,
         "height": 90,
         "x": 24,
         "y": 51
        },
        {
         "id": "B",
         "labels": [
          {
           "text": "B",
           "width": 10,
           "height": 12,
           "x": 155,
           "y": 24
          }
         ],
         "ports": [
          {
           "id": "p",
           "width": 10,
           "height": 10,
           "labels": [
            {
             "text": "p",
             "width": 10,
             "height": 12,
             "x": 134,
             "y": 63
            }
           ],
           "x": 145,
           "y": 76
          }
         ],
         "children": [
          {
           "id": "y",
           "width": 50,
           "height": 90,
           "x": 167,
           "y": 36
          },
          {
           "id": "z",
           "width": 50,
           "height": 90,
           "x": 237,
           "y": 51
          }
         ],
         "edges": [
          {
           "id": "e1",
           "sources": [
            "y"
           ],
           "targets": [
            "z"
           ],
           "sections": [
            {
             "id": "e1_s0",
             "startPoint": {
              "x": 62,
              "y": 57
             },
             "endPoint": {
              "x": 82,
              "y": 57
             },
             "incomingShape": "y",
             "outgoingShape": "z"
            }
           ],
           "container": "B"
          },
          {
           "id": "e2",
           "sources": [
            "x"
           ],
           "targets": [
            "z"
           ],
           "labels": [
            {
             "text": "e2",
             "width": 20,
             "height": 12,
             "x": 82,
             "y": 102
            }
           ],
           "sections": [
            {
             "id": "e2_s0",
             "startPoint": {
              "x": 62,
              "y": 99
             },
             "endPoint": {
              "x": 225,
              "y": 99
             },
             "bendPoints": [
              {
               "x": 112,
               "y": 99
              },
              {
               "x": 112,
               "y": 124
              },
              {
               "x": 215,
               "y": 124
              },
              {
               "x": 215,
               "y": 99
              }
             ],
             "incomingShape": "x",
             "outgoingShape": "z"
            }
           ],
           "container": "A"
          },
          {
           "id": "e3",
           "sources": [
            "x"
           ],
           "targets": [
            "p"
           ],
           "sections": [
            {
             "id": "e3_s0",
             "startPoint": {
              "x": 62,
              "y": 69
             },
             "endPoint": {
              "x": 133,
              "y": 69
             },
             "incomingShape": "x",
             "outgoingShape": "p"
            }
           ],
           "container": "A"
          },
          {
           "id": "e4",
           "sources": [
            "p"
           ],
           "targets": [
            "y"
           ],
           "sections": [
            {
             "id": "e4_s0",
             "startPoint": {
              "x": -10,
              "y": 57
             },
             "endPoint": {
              "x": 12,
              "y": 57
             },
             "incomingShape": "p",
             "outgoingShape": "y"
            }
           ],
           "container": "B"
          }
         ],
         "x": 155,
         "y": 24,
         "width": 144,
         "height": 129
        }
       ],
       "x": 12,
       "y": 12,
       "width": 299,
       "height": 153
      }
     ],
     "x": 0,
     "y": 0,
     "width": 323,
     "height": 177
    }
    '''
    
    val outputRP = '''
    {
     "id": "root",
     "properties": {
      "algorithm": "layered",
      "org.eclipse.elk.hierarchyHandling": "INCLUDE_CHILDREN",
      "org.eclipse.elk.json.shapeCoords": "ROOT",
      "org.eclipse.elk.json.edgeCoords": "PARENT"
     },
     "children": [
      {
       "id": "A",
       "children": [
        {
         "id": "x",
         "width": 50,
         "height": 90,
         "x": 24,
         "y": 51
        },
        {
         "id": "B",
         "labels": [
          {
           "text": "B",
           "width": 10,
           "height": 12,
           "x": 155,
           "y": 24
          }
         ],
         "ports": [
          {
           "id": "p",
           "width": 10,
           "height": 10,
           "labels": [
            {
             "text": "p",
             "width": 10,
             "height": 12,
             "x": 134,
             "y": 63
            }
           ],
           "x": 145,
           "y": 76
          }
         ],
         "children": [
          {
           "id": "y",
           "width": 50,
           "height": 90,
           "x": 167,
           "y": 36
          },
          {
           "id": "z",
           "width": 50,
           "height": 90,
           "x": 237,
           "y": 51
          }
         ],
         "edges": [
          {
           "id": "e1",
           "sources": [
            "y"
           ],
           "targets": [
            "z"
           ],
           "sections": [
            {
             "id": "e1_s0",
             "startPoint": {
              "x": 62,
              "y": 57
             },
             "endPoint": {
              "x": 82,
              "y": 57
             },
             "incomingShape": "y",
             "outgoingShape": "z"
            }
           ]
          },
          {
           "id": "e2",
           "sources": [
            "x"
           ],
           "targets": [
            "z"
           ],
           "labels": [
            {
             "text": "e2",
             "width": 20,
             "height": 12,
             "x": -61,
             "y": 90
            }
           ],
           "sections": [
            {
             "id": "e2_s0",
             "startPoint": {
              "x": -81,
              "y": 87
             },
             "endPoint": {
              "x": 82,
              "y": 87
             },
             "bendPoints": [
              {
               "x": -31,
               "y": 87
              },
              {
               "x": -31,
               "y": 112
              },
              {
               "x": 72,
               "y": 112
              },
              {
               "x": 72,
               "y": 87
              }
             ],
             "incomingShape": "x",
             "outgoingShape": "z"
            }
           ]
          },
          {
           "id": "e3",
           "sources": [
            "x"
           ],
           "targets": [
            "p"
           ],
           "sections": [
            {
             "id": "e3_s0",
             "startPoint": {
              "x": -81,
              "y": 57
             },
             "endPoint": {
              "x": -10,
              "y": 57
             },
             "incomingShape": "x",
             "outgoingShape": "p"
            }
           ]
          },
          {
           "id": "e4",
           "sources": [
            "p"
           ],
           "targets": [
            "y"
           ],
           "sections": [
            {
             "id": "e4_s0",
             "startPoint": {
              "x": -10,
              "y": 57
             },
             "endPoint": {
              "x": 12,
              "y": 57
             },
             "incomingShape": "p",
             "outgoingShape": "y"
            }
           ]
          }
         ],
         "x": 155,
         "y": 24,
         "width": 144,
         "height": 129
        }
       ],
       "x": 12,
       "y": 12,
       "width": 299,
       "height": 153
      }
     ],
     "x": 0,
     "y": 0,
     "width": 323,
     "height": 177
    }
    
    '''
    
    val outputRR = '''
    {
     "id": "root",
     "properties": {
      "algorithm": "layered",
      "org.eclipse.elk.hierarchyHandling": "INCLUDE_CHILDREN",
      "org.eclipse.elk.json.shapeCoords": "ROOT",
      "org.eclipse.elk.json.edgeCoords": "ROOT"
     },
     "children": [
      {
       "id": "A",
       "children": [
        {
         "id": "x",
         "width": 50,
         "height": 90,
         "x": 24,
         "y": 51
        },
        {
         "id": "B",
         "labels": [
          {
           "text": "B",
           "width": 10,
           "height": 12,
           "x": 155,
           "y": 24
          }
         ],
         "ports": [
          {
           "id": "p",
           "width": 10,
           "height": 10,
           "labels": [
            {
             "text": "p",
             "width": 10,
             "height": 12,
             "x": 134,
             "y": 63
            }
           ],
           "x": 145,
           "y": 76
          }
         ],
         "children": [
          {
           "id": "y",
           "width": 50,
           "height": 90,
           "x": 167,
           "y": 36
          },
          {
           "id": "z",
           "width": 50,
           "height": 90,
           "x": 237,
           "y": 51
          }
         ],
         "edges": [
          {
           "id": "e1",
           "sources": [
            "y"
           ],
           "targets": [
            "z"
           ],
           "sections": [
            {
             "id": "e1_s0",
             "startPoint": {
              "x": 217,
              "y": 81
             },
             "endPoint": {
              "x": 237,
              "y": 81
             },
             "incomingShape": "y",
             "outgoingShape": "z"
            }
           ]
          },
          {
           "id": "e2",
           "sources": [
            "x"
           ],
           "targets": [
            "z"
           ],
           "labels": [
            {
             "text": "e2",
             "width": 20,
             "height": 12,
             "x": 94,
             "y": 114
            }
           ],
           "sections": [
            {
             "id": "e2_s0",
             "startPoint": {
              "x": 74,
              "y": 111
             },
             "endPoint": {
              "x": 237,
              "y": 111
             },
             "bendPoints": [
              {
               "x": 124,
               "y": 111
              },
              {
               "x": 124,
               "y": 136
              },
              {
               "x": 227,
               "y": 136
              },
              {
               "x": 227,
               "y": 111
              }
             ],
             "incomingShape": "x",
             "outgoingShape": "z"
            }
           ]
          },
          {
           "id": "e3",
           "sources": [
            "x"
           ],
           "targets": [
            "p"
           ],
           "sections": [
            {
             "id": "e3_s0",
             "startPoint": {
              "x": 74,
              "y": 81
             },
             "endPoint": {
              "x": 145,
              "y": 81
             },
             "incomingShape": "x",
             "outgoingShape": "p"
            }
           ]
          },
          {
           "id": "e4",
           "sources": [
            "p"
           ],
           "targets": [
            "y"
           ],
           "sections": [
            {
             "id": "e4_s0",
             "startPoint": {
              "x": 145,
              "y": 81
             },
             "endPoint": {
              "x": 167,
              "y": 81
             },
             "incomingShape": "p",
             "outgoingShape": "y"
            }
           ]
          }
         ],
         "x": 155,
         "y": 24,
         "width": 144,
         "height": 129
        }
       ],
       "x": 12,
       "y": 12,
       "width": 299,
       "height": 153
      }
     ],
     "x": 0,
     "y": 0,
     "width": 323,
     "height": 177
    }
    '''
}
