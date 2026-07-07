/**
 * ******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *  ******************************************************************************
 */
package org.eclipse.elk.core.debug.grandom.gRandom;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.eclipse.emf.common.util.Enumerator;

/**
 * <!-- begin-user-doc -->
 * A representation of the literals of the enumeration '<em><b>Form</b></em>',
 * and utility methods for working with them.
 * <!-- end-user-doc -->
 * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getForm()
 * @model
 * @generated
 */
public enum Form implements Enumerator
{
  /**
   * The '<em><b>Trees</b></em>' literal object.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #TREES_VALUE
   * @generated
   * @ordered
   */
  TREES(0, "trees", "trees"),

  /**
   * The '<em><b>Custom</b></em>' literal object.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #CUSTOM_VALUE
   * @generated
   * @ordered
   */
  CUSTOM(1, "custom", "graphs"),

  /**
   * The '<em><b>Bipartite</b></em>' literal object.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #BIPARTITE_VALUE
   * @generated
   * @ordered
   */
  BIPARTITE(2, "bipartite", "bipartite graphs"),

  /**
   * The '<em><b>Biconnected</b></em>' literal object.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #BICONNECTED_VALUE
   * @generated
   * @ordered
   */
  BICONNECTED(3, "biconnected", "biconnected graphs"),

  /**
   * The '<em><b>Triconnected</b></em>' literal object.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #TRICONNECTED_VALUE
   * @generated
   * @ordered
   */
  TRICONNECTED(4, "triconnected", "triconnected graphs"),

  /**
   * The '<em><b>Acyclic</b></em>' literal object.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #ACYCLIC_VALUE
   * @generated
   * @ordered
   */
  ACYCLIC(5, "acyclic", "acyclic graphs");

  /**
   * The '<em><b>Trees</b></em>' literal value.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of '<em><b>Trees</b></em>' literal object isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @see #TREES
   * @model name="trees"
   * @generated
   * @ordered
   */
  public static final int TREES_VALUE = 0;

  /**
   * The '<em><b>Custom</b></em>' literal value.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of '<em><b>Custom</b></em>' literal object isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @see #CUSTOM
   * @model name="custom" literal="graphs"
   * @generated
   * @ordered
   */
  public static final int CUSTOM_VALUE = 1;

  /**
   * The '<em><b>Bipartite</b></em>' literal value.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of '<em><b>Bipartite</b></em>' literal object isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @see #BIPARTITE
   * @model name="bipartite" literal="bipartite graphs"
   * @generated
   * @ordered
   */
  public static final int BIPARTITE_VALUE = 2;

  /**
   * The '<em><b>Biconnected</b></em>' literal value.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of '<em><b>Biconnected</b></em>' literal object isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @see #BICONNECTED
   * @model name="biconnected" literal="biconnected graphs"
   * @generated
   * @ordered
   */
  public static final int BICONNECTED_VALUE = 3;

  /**
   * The '<em><b>Triconnected</b></em>' literal value.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of '<em><b>Triconnected</b></em>' literal object isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @see #TRICONNECTED
   * @model name="triconnected" literal="triconnected graphs"
   * @generated
   * @ordered
   */
  public static final int TRICONNECTED_VALUE = 4;

  /**
   * The '<em><b>Acyclic</b></em>' literal value.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of '<em><b>Acyclic</b></em>' literal object isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @see #ACYCLIC
   * @model name="acyclic" literal="acyclic graphs"
   * @generated
   * @ordered
   */
  public static final int ACYCLIC_VALUE = 5;

  /**
   * An array of all the '<em><b>Form</b></em>' enumerators.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private static final Form[] VALUES_ARRAY =
    new Form[]
    {
      TREES,
      CUSTOM,
      BIPARTITE,
      BICONNECTED,
      TRICONNECTED,
      ACYCLIC,
    };

  /**
   * A public read-only list of all the '<em><b>Form</b></em>' enumerators.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public static final List<Form> VALUES = Collections.unmodifiableList(Arrays.asList(VALUES_ARRAY));

  /**
   * Returns the '<em><b>Form</b></em>' literal with the specified literal value.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param literal the literal.
   * @return the matching enumerator or <code>null</code>.
   * @generated
   */
  public static Form get(String literal)
  {
    for (int i = 0; i < VALUES_ARRAY.length; ++i)
    {
      Form result = VALUES_ARRAY[i];
      if (result.toString().equals(literal))
      {
        return result;
      }
    }
    return null;
  }

  /**
   * Returns the '<em><b>Form</b></em>' literal with the specified name.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param name the name.
   * @return the matching enumerator or <code>null</code>.
   * @generated
   */
  public static Form getByName(String name)
  {
    for (int i = 0; i < VALUES_ARRAY.length; ++i)
    {
      Form result = VALUES_ARRAY[i];
      if (result.getName().equals(name))
      {
        return result;
      }
    }
    return null;
  }

  /**
   * Returns the '<em><b>Form</b></em>' literal with the specified integer value.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the integer value.
   * @return the matching enumerator or <code>null</code>.
   * @generated
   */
  public static Form get(int value)
  {
    switch (value)
    {
      case TREES_VALUE: return TREES;
      case CUSTOM_VALUE: return CUSTOM;
      case BIPARTITE_VALUE: return BIPARTITE;
      case BICONNECTED_VALUE: return BICONNECTED;
      case TRICONNECTED_VALUE: return TRICONNECTED;
      case ACYCLIC_VALUE: return ACYCLIC;
    }
    return null;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private final int value;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private final String name;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private final String literal;

  /**
   * Only this class can construct instances.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private Form(int value, String name, String literal)
  {
    this.value = value;
    this.name = name;
    this.literal = literal;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public int getValue()
  {
    return value;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public String getName()
  {
    return name;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public String getLiteral()
  {
    return literal;
  }

  /**
   * Returns the literal value of the enumerator, which is its string representation.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public String toString()
  {
    return literal;
  }
  
} //Form
